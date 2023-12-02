// header for ADMM

#pragma once
// STL
#include <complex>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

// Additional libraries
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Structures
namespace ADMM{
    static inline double pi(){
        double value = 3.14159265358979323846;
        return value;
    }

    struct opf_struct{
        // System statistic
        struct statistic_struct{
            int num_node;
            int num_line;
            int num_variable;
            int num_constraint;
        };
        statistic_struct statistic;

        // Network information
        struct network_struct{
            // Network topology information
            std::vector <Eigen::Vector2i> topology;

            // Line conductance information
            Eigen::VectorXcd line_conductance;
        };
        network_struct network;

        // Objective function
        struct obj_struct{
            struct cost_func_struct{
                // Structures for bids
                struct bid_struct{
                    Eigen::VectorXd price;
                    Eigen::VectorXd quantity;
                };
                bid_struct demand;
                bid_struct supply;

                // Substructure for merit order curve
                struct moc_struct{
                    Eigen::VectorXd price;
                    Eigen::VectorXd quantity;
                    Eigen::VectorXd obj;
                };
                moc_struct moc;

                void moc_set(){
                    int num_row = 2 * (this->demand.price.size() + this->supply.price.size() - 3);
                    this->moc.price = Eigen::VectorXd(num_row);
                    this->moc.quantity = Eigen::VectorXd(num_row);
                    this->moc.obj = Eigen::VectorXd(num_row);

                    // Initialization for the loop
                    int price_demand_ID = 0;
                    int price_supply_ID = 0;
                    int price_ID = 0;
                    double obj_value = -(this->demand.price.segment(1, this->demand.price.size() - 2).array() * this->demand.quantity.segment(1, this->demand.quantity.size() - 2).array()).sum();

                    // First item when price = -Inf
                    this->moc.price(0) = -std::numeric_limits<double>::infinity();
                    this->moc.quantity(0) = -this->demand.quantity.sum();
                    this->moc.obj(0) = obj_value;

                    while(price_ID < num_row - 2){
                        price_ID += 1;

                        // Next demand marginal price lower than supply; add that demand quantity to residual load
                        if(this->demand.price(price_demand_ID + 1) <= this->supply.price(price_supply_ID + 1)){
                            price_demand_ID += price_demand_ID < this->demand.price.size() - 2;
                            this->moc.price(price_ID) = this->demand.price(price_demand_ID);
                            this->moc.quantity(price_ID) = this->moc.quantity(price_ID - 1);
                            this->moc.obj(price_ID, 2) = obj_value;

                            price_ID += 1;
                            obj_value += this->demand.price(price_demand_ID) * this->demand.quantity(price_demand_ID);
                            this->moc.price(price_ID) = this->demand.price(price_demand_ID);
                            this->moc.quantity(price_ID) = this->moc.quantity(price_ID - 1) + this->demand.quantity(price_demand_ID);
                            this->moc.obj(price_ID) = obj_value;
                        }
                        // Next supply marginal price lower than demand; add that supply quantity to residual load
                        else{
                            price_supply_ID += price_supply_ID < this->supply.price.size() - 2;
                            this->moc.price(price_ID) = this->supply.price(price_supply_ID);
                            this->moc.quantity(price_ID) = this->moc.quantity(price_ID - 1);
                            this->moc.obj(price_ID, 2) = obj_value;

                            price_ID += 1;
                            obj_value += this->supply.price(price_supply_ID) * this->supply.quantity(price_supply_ID);
                            this->moc.price(price_ID) = this->supply.price(price_supply_ID);
                            this->moc.quantity(price_ID) = this->moc.quantity(price_ID - 1) + this->supply.quantity(price_supply_ID);
                            this->moc.obj(price_ID) = obj_value;
                        }
                    }

                    // Last item when price = Inf
                    price_ID += 1;
                    this->moc.price(price_ID) = std::numeric_limits<double>::infinity();
                    this->moc.quantity(price_ID) = this->moc.quantity(price_ID - 1);
                    this->moc.obj(price_ID) = obj_value;
                }
            };

            struct transform_struct{
                Eigen::VectorXd scale;
                Eigen::VectorXd shift;
                double obj;
            };
            transform_struct transformation;

            std::vector <cost_func_struct> cost_funcs;
        };
        obj_struct obj;

        // Solver
        struct solver_struct{
            // Augmented Lagrangian:
            // L(x, u) = f(x) + .5 * \rho * ||A %*% x - c + u||^2
            // x = scaled voltage (V), power source / sink (S), and line current (I)
            // c = equality boundary
            // u = dual variables of the constraints
            // f(x): cost function associated with V, S, and I (f(V) and f(I) are well functions)

            struct constraint_struct{
                // Main matrix
                Eigen::SparseMatrix <double> Matrix_main;
                Eigen::SparseMatrix <double> PSD_main;
                std::vector <std::vector <std::pair <int, double>>> Mat_2_main_terms;
                std::vector <std::vector <std::pair <int, double>>> PSD_main_terms;
                Eigen::VectorXd sensitivity_var;
                Eigen::VectorXd boundary;
            };
            constraint_struct constraint;

            // Solution
            struct sol_struct{
                struct value_struct{
                    struct variable_struct{
                        Eigen::VectorXd now;
                        Eigen::VectorXd prev;
                    };
                    variable_struct variables;
                    Eigen::VectorXd price_margin;
                    Eigen::VectorXd error;
                };

                value_struct prime;
                value_struct dual;
                double obj_value;
            };
            sol_struct sol;
        };
        solver_struct solver;

        // Functions
        void problem_size_parmeters_set(int factor = 1){
            this->statistic.num_variable = factor * (2 * this->statistic.num_node + this->statistic.num_line);
            this->statistic.num_constraint = factor * (this->statistic.num_node + this->statistic.num_line);
        }

        void transformation_set(bool flag = 1){
            // transformation:
            // x' = 1 / m * (x - x_0)
            // price and quantity in moc must be scaled accordingly
            // price: p(x') = m * p(x)
            this->obj.transformation.scale = Eigen::VectorXd::Ones(this->statistic.num_variable);
            this->obj.transformation.shift = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->obj.transformation.obj = 1.;
            if(!flag){
                return;
            }

            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->obj.transformation.scale(var_iter) = this->obj.cost_funcs[var_iter].moc.quantity.maxCoeff();
                this->obj.transformation.scale(var_iter) -= this->obj.cost_funcs[var_iter].moc.quantity.minCoeff();
                this->obj.transformation.scale(var_iter) /= 2;
                this->obj.transformation.shift(var_iter) = this->obj.cost_funcs[var_iter].moc.quantity.maxCoeff();
                this->obj.transformation.shift(var_iter) += this->obj.cost_funcs[var_iter].moc.quantity.minCoeff();
                this->obj.transformation.shift(var_iter) /= 2;

                this->obj.cost_funcs[var_iter].moc.quantity -= this->obj.transformation.shift(var_iter) * Eigen::VectorXd::Ones(this->statistic.num_variable);
                this->obj.cost_funcs[var_iter].moc.quantity /= this->obj.transformation.scale(var_iter);
                this->obj.cost_funcs[var_iter].moc.price *= this->obj.transformation.scale(var_iter);
            }

            this->obj.transformation.obj = 0.;
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->obj.transformation.obj = std::max(this->obj.transformation.obj, abs(this->obj.cost_funcs[var_iter].moc.price[this->obj.cost_funcs[var_iter].moc.price.size() - 2]));
            }
            this->obj.transformation.obj *= 1E2;
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->obj.cost_funcs[var_iter].moc.price = this->obj.cost_funcs[var_iter].moc.price / this->obj.transformation.obj;
            }
        }

        // Main matrix initialization for DC OPF
        void DC_Matrix_main_set(){
            // Variables: V, S, I
            // Constraints: Node balance, line current
            this->solver.constraint.Mat_2_main_terms = std::vector <std::vector <std::pair <int, double>>> (this->statistic.num_variable);
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->solver.constraint.Mat_2_main_terms[var_iter].reserve(this->statistic.num_constraint);
            }
            this->solver.constraint.PSD_main_terms = std::vector <std::vector <std::pair <int, double>>> (this->statistic.num_variable);
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->solver.constraint.PSD_main_terms[var_iter].reserve(this->statistic.num_variable);
            }
            this->solver.constraint.sensitivity_var = Eigen::VectorXd (this->statistic.num_variable);

            // Set the main matrix
            Eigen::MatrixXd Mat = Eigen::MatrixXd::Zero(this->statistic.num_constraint, this->statistic.num_variable);

            // Node Balance Equation
            // S - t(NL) %*% I = 0
            for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                int var_ID = this->statistic.num_node + node_iter;
                Mat(node_iter, var_ID) = 1.;
            }
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = 2 * this->statistic.num_node + line_iter;
                Mat(this->network.topology[line_iter](0), var_ID) = -1.;
                Mat(this->network.topology[line_iter](1), var_ID) = 1.;
            }

            // Line Current Equation
            // Y_l %*% NL %*% V - I = 0
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = 2 * this->statistic.num_node + line_iter;
                int constr_ID = this->statistic.num_node + line_iter;
                double y_l = this->network.line_conductance(line_iter).imag();

                Mat(constr_ID, this->network.topology[line_iter](0)) = y_l;
                Mat(constr_ID, this->network.topology[line_iter](1)) = -y_l;
                Mat(constr_ID, var_ID) = -1.;
            }

            // Set boundary for equality constraints
            this->solver.constraint.boundary = -Mat * this->obj.transformation.shift;

            // Apply transformation to the matrix
            Eigen::MatrixXd Diag = Eigen::MatrixXd::Zero(this->statistic.num_variable, this->statistic.num_variable);
            Eigen::MatrixXd Diag_inv = Eigen::MatrixXd::Zero(this->statistic.num_variable, this->statistic.num_variable);
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                Diag(var_iter, var_iter) = this->obj.transformation.scale(var_iter);
                Diag_inv(var_iter, var_iter) = 1. / Diag(var_iter, var_iter);
            }
            Mat = Mat * Diag;
            Eigen::MatrixXd PSD = Mat.transpose() * Mat;

            // Record non-zero terms for Mat
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                for(int constr_iter = 0; constr_iter < this->statistic.num_constraint; ++ constr_iter){
                    if(Mat(constr_iter, var_iter) != 0.){
                        this->solver.constraint.Mat_2_main_terms[var_iter].push_back(std::pair <int, double> (constr_iter, Mat(constr_iter, var_iter)));
                    }
                }
            }

            // Record non-zero terms for PSD
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                for(int row_iter = 0; row_iter < this->statistic.num_variable; ++ row_iter){
                    if(PSD(row_iter, var_iter) != 0.){
                        if(row_iter == var_iter){
                            this->solver.constraint.sensitivity_var(var_iter) = PSD(row_iter, var_iter);
                        }
                        else{
                            this->solver.constraint.PSD_main_terms[var_iter].push_back(std::pair <int, double> (row_iter, PSD(row_iter, var_iter)));
                        }
                    }
                }
            }

            ////////////////////////////////////////////////////////////////////////
            int num_entry = this->statistic.num_node + 5 * this->statistic.num_line;

            // Set the main matrix
            this->solver.constraint.Matrix_main = Eigen::SparseMatrix <double> (this->statistic.num_constraint, this->statistic.num_variable);
            std::vector <Eigen::Triplet <double>> Matrix_main_trip;
            Matrix_main_trip.reserve(num_entry);

            // Node Balance Equation
            // S - t(NL) %*% I = 0
            for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                int var_ID = this->statistic.num_node + node_iter;
                Matrix_main_trip.push_back(Eigen::Triplet <double> (node_iter, var_ID, 1));
            }
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = 2 * this->statistic.num_node + line_iter;
                Matrix_main_trip.push_back(Eigen::Triplet <double> (this->network.topology[line_iter](0), var_ID, -1));
                Matrix_main_trip.push_back(Eigen::Triplet <double> (this->network.topology[line_iter](1), var_ID, 1));
            }

            // Line Current Equation
            // Y_l %*% NL %*% V - I = 0
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = 2 * this->statistic.num_node + line_iter;
                int constr_ID = this->statistic.num_node + line_iter;
                double y_l = this->network.line_conductance(line_iter).imag();

                Matrix_main_trip.push_back(Eigen::Triplet <double> (constr_ID, this->network.topology[line_iter](0), y_l));
                Matrix_main_trip.push_back(Eigen::Triplet <double> (constr_ID, this->network.topology[line_iter](1), -y_l));
                Matrix_main_trip.push_back(Eigen::Triplet <double> (constr_ID, var_ID, -1));
            }

            // Put Everything in to the matrix
            this->solver.constraint.Matrix_main.setFromTriplets(Matrix_main_trip.begin(), Matrix_main_trip.end());

            // Set boundary for equality constraints
            this->solver.constraint.boundary = -this->solver.constraint.Matrix_main * this->obj.transformation.shift;

            // Apply transformation to the matrix
            Eigen::SparseMatrix <double> Diagonal(this->statistic.num_variable, this->statistic.num_variable);
            std::vector <Eigen::Triplet <double>> Diagonal_trip;
            Diagonal_trip.reserve(this->statistic.num_variable);
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                Diagonal_trip.push_back(Eigen::Triplet <double> (var_iter, var_iter, this->obj.transformation.scale(var_iter)));
            }
            Diagonal.setFromTriplets(Diagonal_trip.begin(), Diagonal_trip.end());
            this->solver.constraint.Matrix_main = this->solver.constraint.Matrix_main * Diagonal;
            this->solver.constraint.PSD_main = this->solver.constraint.Matrix_main.transpose() * this->solver.constraint.Matrix_main;
        }

        // Solver function: direct search the root
        void solve_iteration_subproblem(double rho){
            // Update S and V
            for(int var_iter = 0; var_iter < 2 * this->statistic.num_node; ++ var_iter){
                double inner_prod = 0.;
                // add from PSD
                int entry_iter = 0;
                while(entry_iter < this->solver.constraint.PSD_main_terms[var_iter].size()){
                    int row_ID = this->solver.constraint.PSD_main_terms[var_iter][entry_iter].first;
                    double coeff = this->solver.constraint.PSD_main_terms[var_iter][entry_iter].second;
                    double value = coeff * this->solver.sol.prime.variables.now(row_ID);
                    inner_prod += value;
                    entry_iter += 1;
                }

                // add from Mat
                entry_iter = 0;
                while(entry_iter < this->solver.constraint.Mat_2_main_terms[var_iter].size()){
                    int constr_ID = this->solver.constraint.Mat_2_main_terms[var_iter][entry_iter].first;
                    double coeff = this->solver.constraint.Mat_2_main_terms[var_iter][entry_iter].second;
                    double value = coeff * (this->solver.sol.dual.variables.now(constr_ID) - this->solver.constraint.boundary(constr_ID));
                    inner_prod += value;
                    entry_iter += 1;
                }

                Eigen::Vector2i gap_price_ID(0, this->obj.cost_funcs[var_iter].moc.price.size() - 1);
                int mid_price_ID = gap_price_ID.sum() / 2;
                while(gap_price_ID(1) - gap_price_ID(0) > 1){
                    Eigen::Vector2d gap_rent;
                    gap_rent(0) = this->obj.cost_funcs[var_iter].moc.price(gap_price_ID(0));
                    gap_rent(0) += rho * (this->solver.constraint.sensitivity_var(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(0)) + inner_prod);
                    gap_rent(1) = this->obj.cost_funcs[var_iter].moc.price(gap_price_ID(1));
                    gap_rent(1) += rho * (this->solver.constraint.sensitivity_var(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(1)) + inner_prod);

                    double mid_rent;
                    mid_rent = this->obj.cost_funcs[var_iter].moc.price(mid_price_ID);
                    mid_rent += rho * (this->solver.constraint.sensitivity_var(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID) + inner_prod);

                    // Bisection root-search method
                    if(mid_rent * gap_rent(0) <= 0){
                        gap_price_ID(1) = mid_price_ID;
                    }
                    else{
                        gap_price_ID(0) = mid_price_ID;
                    }

                    mid_price_ID = gap_price_ID.sum() / 2;
                }

                // Check binding term
                if(gap_price_ID(0) % 2 == 1){
                    // price fixed case
                    this->solver.sol.prime.price_margin(var_iter) = this->obj.cost_funcs[var_iter].moc.price(mid_price_ID);
                    this->solver.sol.prime.variables.now(var_iter) = -this->solver.sol.prime.price_margin(var_iter) / rho;
                    this->solver.sol.prime.variables.now(var_iter) -= inner_prod;
                    this->solver.sol.prime.variables.now(var_iter) /= this->solver.constraint.sensitivity_var(var_iter);
                }
                else{
                    // quantity fixed case
                    this->solver.sol.prime.variables.now(var_iter) = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID);
                    this->solver.sol.prime.price_margin(var_iter) = -rho * (this->solver.constraint.sensitivity_var(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID) + inner_prod);
                }
            }

            // Update I
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int I_ID = 2 * this->statistic.num_node + line_iter;

                double inner_prod = 0.;
                // add from PSD
                int entry_iter = 0;
                while(entry_iter < this->solver.constraint.PSD_main_terms[I_ID].size()){
                    int row_ID = this->solver.constraint.PSD_main_terms[I_ID][entry_iter].first;
                    double coeff = this->solver.constraint.PSD_main_terms[I_ID][entry_iter].second;
                    double value = coeff * this->solver.sol.prime.variables.now(row_ID);;
                    inner_prod += value;
                    entry_iter += 1;
                }

                // add from Mat
                entry_iter = 0;
                while(entry_iter < this->solver.constraint.Mat_2_main_terms[I_ID].size()){
                    int constr_ID = this->solver.constraint.Mat_2_main_terms[I_ID][entry_iter].first;
                    double coeff = this->solver.constraint.Mat_2_main_terms[I_ID][entry_iter].second;
                    double value = coeff * (this->solver.sol.dual.variables.now(constr_ID) - this->solver.constraint.boundary(constr_ID));
                    inner_prod += value;
                    entry_iter += 1;
                }

                Eigen::Vector2i gap_price_ID(0, this->obj.cost_funcs[I_ID].moc.price.size() - 1);
                int mid_price_ID = gap_price_ID.sum() / 2;
                while(gap_price_ID(1) - gap_price_ID(0) > 1){
                    Eigen::Vector2d gap_rent;
                    gap_rent(0) = this->obj.cost_funcs[I_ID].moc.price(gap_price_ID(0));
                    gap_rent(0) += rho * (this->solver.constraint.sensitivity_var(I_ID) * this->obj.cost_funcs[I_ID].moc.quantity(gap_price_ID(0)) + inner_prod);
                    gap_rent(1) = this->obj.cost_funcs[I_ID].moc.price(gap_price_ID(1));
                    gap_rent(1) += rho * (this->solver.constraint.sensitivity_var(I_ID) * this->obj.cost_funcs[I_ID].moc.quantity(gap_price_ID(1)) + inner_prod);

                    double mid_rent;
                    mid_rent = this->obj.cost_funcs[I_ID].moc.price(mid_price_ID);
                    mid_rent += rho * (this->solver.constraint.sensitivity_var(I_ID) * this->obj.cost_funcs[I_ID].moc.quantity(mid_price_ID) + inner_prod);

                    // Bisection root-search method
                    if(mid_rent * gap_rent(0) <= 0){
                        gap_price_ID(1) = mid_price_ID;
                    }
                    else{
                        gap_price_ID(0) = mid_price_ID;
                    }

                    mid_price_ID = gap_price_ID.sum() / 2;
                }

                // Check binding term
                if(gap_price_ID(0) % 2 == 1){
                    // price fixed case
                    this->solver.sol.prime.price_margin(I_ID) = this->obj.cost_funcs[I_ID].moc.price(mid_price_ID);
                    this->solver.sol.prime.variables.now(I_ID) = -this->solver.sol.prime.price_margin(I_ID) / rho;
                    this->solver.sol.prime.variables.now(I_ID) -= inner_prod;
                    this->solver.sol.prime.variables.now(I_ID) /= this->solver.constraint.sensitivity_var(I_ID);
                }
                else{
                    // quantity fixed case
                    this->solver.sol.prime.variables.now(I_ID) = this->obj.cost_funcs[I_ID].moc.quantity(mid_price_ID);
                    this->solver.sol.prime.price_margin(I_ID) = -rho * (this->solver.constraint.sensitivity_var(I_ID) * this->obj.cost_funcs[I_ID].moc.quantity(mid_price_ID) + inner_prod);
                }
            }
        }

        void solve_root(double tol_prime, double tol_dual, bool print_flag = 1){
            // Initialization
            double rho = 1.;
            double omega = 1.;
            double alpha = 1.;
            this->solver.sol.prime.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.variables.now = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.price_margin = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.dual.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_constraint);
            this->solver.sol.dual.variables.now = Eigen::VectorXd::Zero(this->statistic.num_constraint);

            // Main loop
            // f(x) + .5 * rho * ||A %*% x - c + u||^2
            // f(x) + .5 * rho * (M_{ij} x_i x_j + 2 * A_{ij} b_i x_j + ...)
            // KKT: p_i + rho * (M_{ij} x_j + t(A)_{ij} b_j) = 0
            int loop = 0;
            while(true){
                // Update prime variables and associate prices
                int sub_loop = 0;
                while(true){
                    solve_iteration_subproblem(rho);

                    // Stabilization weight for x update
                    this->solver.sol.prime.variables.now = omega * this->solver.sol.prime.variables.now + (1. - omega) * this->solver.sol.prime.variables.prev;
                    this->solver.sol.prime.variables.prev = this->solver.sol.prime.variables.now;

                    // Check convergence for subproblem
                    double tol_sub = exp(-loop);
                    tol_sub = std::max(tol_sub, 1E-6 * tol_dual);

                    this->solver.sol.dual.error = this->solver.constraint.Matrix_main * this->solver.sol.prime.variables.now - this->solver.constraint.boundary + this->solver.sol.dual.variables.now;
                    this->solver.sol.dual.error = rho * this->solver.constraint.Matrix_main.transpose() * this->solver.sol.dual.error;
                    this->solver.sol.dual.error += this->solver.sol.prime.price_margin;
                    if((this->solver.sol.dual.error.array() * this->solver.sol.dual.error.array()).maxCoeff() < tol_sub){
                        break;
                    }

                    sub_loop += 1;
                }
//                std::cout << sub_loop << "\n";

                // Update dual variables
                // u <- u + alpha * (A %*% x - c)
                this->solver.sol.dual.variables.now += alpha * (this->solver.constraint.Matrix_main * this->solver.sol.prime.variables.now - this->solver.constraint.boundary);

                // Check convergence for whole problem
                this->solver.sol.prime.error = this->solver.constraint.Matrix_main * this->solver.sol.prime.variables.now - this->solver.constraint.boundary;
                this->solver.sol.dual.error = this->solver.sol.prime.price_margin;
                this->solver.sol.dual.error += rho * this->solver.constraint.Matrix_main.transpose() * this->solver.sol.dual.variables.now;
                double prime_error_norm = this->solver.sol.prime.error.norm();
                prime_error_norm /= this->solver.sol.prime.error.size();
                double dual_error_norm = this->solver.sol.dual.error.norm();
                dual_error_norm /= this->solver.sol.dual.error.size();
                if(prime_error_norm < tol_prime && dual_error_norm < tol_dual){
                    break;
                }

                // Print progress
                if(print_flag){
                    int sys_show = (int) 7. - log(this->statistic.num_variable) / log(10.);
                    sys_show = std::max(sys_show, 0);
//                    if(loop % (int) pow(10., sys_show) == 0){
                    if(loop % 1 == 0){
                        std::cout << "Loop:\t" << loop << "\n";
                        std::cout << "Prime Error:\t" << prime_error_norm << "\n";
                        std::cout << "Dual Error:\t" << dual_error_norm << "\n";
                        std::cout << "\n\n";
                    }
                }
                loop += 1;
            }

            if(print_flag){
                std::cout << "Total loop:\t" << loop << "\n";
                std::cout << "Prime Error:\t" << this->solver.sol.prime.error.norm() / this->solver.sol.prime.error.size() << "\n";
                std::cout << "Dual Error:\t" << this->solver.sol.dual.error.norm() / this->solver.sol.dual.error.size() << "\n";
                std::cout << "Solution:\n";
                std::cout << (this->obj.transformation.scale.array() * this->solver.sol.prime.variables.now.array() + this->obj.transformation.shift.array()).segment(this->statistic.num_node, this->statistic.num_node).transpose() << "\n";
           }
        }

        // solver: gradient method using RK4
        void gradient_x(Eigen::VectorXd &x, Eigen::VectorXd &u, Eigen::VectorXd &dx, double rho){
            // gradient for prime variables
            Eigen::VectorXd price = Eigen::VectorXd::Zero(this->statistic.num_variable);
            for(int var_iter = this->statistic.num_node; var_iter < 2 * this->statistic.num_node; ++ var_iter){
                Eigen::Vector2i gap_price_ID(1, this->obj.cost_funcs[var_iter].moc.price.size() - 2);
                int mid_price_ID = gap_price_ID.sum() / 2;
                while(gap_price_ID(1) - gap_price_ID(0) > 1){
                    Eigen::Vector2d gap_rent;
                    gap_rent(0) = this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(0));
                    gap_rent(0) -= x[var_iter];
                    gap_rent(1) = this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(1));
                    gap_rent(1) -= x[var_iter];

                    double mid_rent;
                    mid_rent = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID);
                    mid_rent -= x[var_iter];

                    // Bisection root-search method
                    if(mid_rent * gap_rent(0) <= 0){
                        gap_price_ID(1) = mid_price_ID;
                    }
                    else{
                        gap_price_ID(0) = mid_price_ID;
                    }

                    mid_price_ID = gap_price_ID.sum() / 2;
                }
                price[var_iter] = this->obj.cost_funcs[var_iter].moc.price(gap_price_ID(0));
            }

            // dx / dt = t(A) * (b - u) - PSD * x
            dx = price / rho + this->solver.constraint.PSD_main * x + this->solver.constraint.Matrix_main.transpose() * (u - this->solver.constraint.boundary);
            dx = -dx;
        }

        void gradient_u(Eigen::VectorXd &x, Eigen::VectorXd &du, double alpha){
            du = this->solver.constraint.Matrix_main * x - this->solver.constraint.boundary;
            du *= alpha;
        }

        void var_constrained(Eigen::VectorXd &x){
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                x[var_iter] = std::max(x[var_iter], this->obj.cost_funcs[var_iter].moc.quantity(0));
                x[var_iter] = std::min(x[var_iter], this->obj.cost_funcs[var_iter].moc.quantity(this->obj.cost_funcs[var_iter].moc.quantity.size() - 1));
            }
        }

        double eigen_value(Eigen::SparseMatrix <double> Mat){
            double eigen_value;
            Eigen::VectorXd eigen_vec_prev = Eigen::VectorXd::Ones(Mat.rows());
            Eigen::VectorXd eigen_vec = eigen_vec_prev;
            while(true){
                eigen_vec = Mat * eigen_vec;
                eigen_value = (eigen_vec.array() * (Mat * eigen_vec).array()).sum();
                eigen_value /= (eigen_vec.array() * eigen_vec.array()).sum();
                eigen_vec /= eigen_value;

                if((eigen_vec - eigen_vec_prev).norm() >= 1E-12){
                    eigen_vec_prev = eigen_vec;
                }
                else{
                    break;
                }
            }

            return(eigen_value);
        }

        void solve_RK4(double tol_prime, double tol_dual, bool print_flag = 1){
            double rho = 1.;
            double dt = 2.7 / eigen_value(this->solver.constraint.PSD_main);
            double alpha = 1. / rho / dt;
            this->solver.sol.prime.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.variables.now = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.price_margin = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.dual.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_constraint);
            this->solver.sol.dual.variables.now = Eigen::VectorXd::Zero(this->statistic.num_constraint);

            int loop = 0;
            while(true){
                Eigen::VectorXd dx_1;
                Eigen::VectorXd du_1;
                gradient_x(this->solver.sol.prime.variables.now, this->solver.sol.dual.variables.now, dx_1, rho);
                gradient_u(this->solver.sol.prime.variables.now, du_1, alpha);

                Eigen::VectorXd x_2 = this->solver.sol.prime.variables.now;
                x_2 += .5 * dx_1 * dt;
                var_constrained(x_2);
                Eigen::VectorXd u_2 = this->solver.sol.dual.variables.now;
                u_2 += .5 * du_1 * dt;
                Eigen::VectorXd dx_2;
                Eigen::VectorXd du_2;
                gradient_x(x_2, u_2, dx_2, rho);
                gradient_u(x_2, du_2, alpha);

                Eigen::VectorXd x_3 = this->solver.sol.prime.variables.now;
                x_3 += .5 * dx_2 * dt;
                var_constrained(x_3);
                Eigen::VectorXd u_3 = this->solver.sol.dual.variables.now;
                u_3 += .5 * du_2 * dt;
                Eigen::VectorXd dx_3;
                Eigen::VectorXd du_3;
                gradient_x(x_3, u_3, dx_3, rho);
                gradient_u(x_3, du_3, alpha);

                Eigen::VectorXd x_4 = this->solver.sol.prime.variables.now;
                x_4 += dx_3 * dt;
                var_constrained(x_4);
                Eigen::VectorXd u_4 = this->solver.sol.dual.variables.now;
                u_4 += du_3 * dt;
                Eigen::VectorXd dx_4;
                gradient_x(x_4, u_4, dx_4, rho);

                this->solver.sol.prime.variables.now += dt / 6 * (dx_1 + 2 * dx_2 + 2 * dx_3 + dx_4);
                var_constrained(this->solver.sol.prime.variables.now);
                Eigen::VectorXd du;
                gradient_u(this->solver.sol.prime.variables.now, du, alpha);
                this->solver.sol.dual.variables.now += du * dt;

                if(loop > 0 && loop % 2000000 == 0){
                    // Root solving algorithm for exact price
                    solve_iteration_subproblem(rho);

                    // Check convergence for whole problem
                    this->solver.sol.prime.error = this->solver.constraint.Matrix_main * this->solver.sol.prime.variables.now - this->solver.constraint.boundary;
                    this->solver.sol.dual.error = this->solver.sol.prime.price_margin;
                    this->solver.sol.dual.error += rho * this->solver.constraint.Matrix_main.transpose() * this->solver.sol.dual.variables.now;
                    double prime_error_norm = this->solver.sol.prime.error.norm();
                    prime_error_norm /= this->solver.sol.prime.error.size();
                    double dual_error_norm = this->solver.sol.dual.error.norm();
                    dual_error_norm /= this->solver.sol.dual.error.size();
                    if(prime_error_norm < tol_prime && dual_error_norm < tol_dual){
                        break;
                    }

                    if(print_flag){
                        int sys_show = (int) 7. - log(this->statistic.num_variable) / log(10.);
                        sys_show = std::max(sys_show, 0);
                        if(loop % (int) pow(10., sys_show) == 0){
                            std::cout << "Loop:\t" << loop << "\n";
                            std::cout << "Prime Error:\t" << prime_error_norm << "\n";
                            std::cout << "Dual Error:\t" << dual_error_norm << "\n";
                            std::cout << "\n\n";
                        }
                    }
                }

                this->solver.sol.prime.variables.prev = this->solver.sol.prime.variables.now;
                this->solver.sol.dual.variables.prev = this->solver.sol.dual.variables.now;
                loop += 1;
            }

            if(print_flag){
                std::cout << "Total loop:\t" << loop << "\n";
                std::cout << "Prime Error:\t" << this->solver.sol.prime.error.norm() / this->solver.sol.prime.error.size() << "\n";
                std::cout << "Dual Error:\t" << this->solver.sol.dual.error.norm() / this->solver.sol.dual.error.size() << "\n";
                std::cout << "Solution:\n";
                std::cout << (this->obj.transformation.scale.array() * this->solver.sol.prime.variables.now.array() + this->obj.transformation.shift.array()).segment(this->statistic.num_node, this->statistic.num_node).transpose() << "\n";
            }
        }

        void solve_RKF45(double tol_prime, double tol_dual, bool print_flag = 1){
            double rho = 1.;
            double dt = 3.1 / eigen_value(this->solver.constraint.PSD_main);
            double alpha = .1 / rho / dt;
            this->solver.sol.prime.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.variables.now = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.price_margin = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.dual.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_constraint);
            this->solver.sol.dual.variables.now = Eigen::VectorXd::Zero(this->statistic.num_constraint);

            int loop = 0;
            while(true){
                Eigen::VectorXd dx_1;
                Eigen::VectorXd du_1;
                gradient_x(this->solver.sol.prime.variables.now, this->solver.sol.dual.variables.now, dx_1, rho);
                gradient_u(this->solver.sol.prime.variables.now, du_1, alpha);

                Eigen::VectorXd x_2 = this->solver.sol.prime.variables.now;
                x_2 += 1. / 4. * dx_1 * dt;
                var_constrained(x_2);
                Eigen::VectorXd u_2 = this->solver.sol.dual.variables.now;
                u_2 += 1. / 4. * du_1 * dt;
                Eigen::VectorXd dx_2;
                Eigen::VectorXd du_2;
                gradient_x(x_2, u_2, dx_2, rho);
                gradient_u(x_2, du_2, alpha);

                Eigen::VectorXd x_3 = this->solver.sol.prime.variables.now;
                x_3 += 3. / 32. * dx_1 * dt + 9. / 32. * dx_2 * dt;
                var_constrained(x_3);
                Eigen::VectorXd u_3 = this->solver.sol.dual.variables.now;
                u_3 += 3. / 32. * du_1 * dt + 9. / 32. * du_2 * dt;
                Eigen::VectorXd dx_3;
                Eigen::VectorXd du_3;
                gradient_x(x_3, u_3, dx_3, rho);
                gradient_u(x_3, du_3, alpha);

                Eigen::VectorXd x_4 = this->solver.sol.prime.variables.now;
                x_4 += 1932. / 2197. * dx_1 * dt - 7200. / 2197. * dx_2 * dt + 7296. / 2197. * dx_3 * dt;
                var_constrained(x_4);
                Eigen::VectorXd u_4 = this->solver.sol.dual.variables.now;
                u_4 += 1932. / 2197. * du_1 * dt - 7200. / 2197. * du_2 * dt + 7296. / 2197. * du_3 * dt;
                Eigen::VectorXd dx_4;
                Eigen::VectorXd du_4;
                gradient_x(x_4, u_4, dx_4, rho);
                gradient_u(x_4, du_4, alpha);

                Eigen::VectorXd x_5 = this->solver.sol.prime.variables.now;
                x_5 += 439. / 216. * dx_1 * dt - 8. * dx_2 * dt + 3680. / 513. * dx_3 * dt - 845. / 4104. * dx_4 * dt;
                var_constrained(x_5);
                Eigen::VectorXd u_5 = this->solver.sol.dual.variables.now;
                u_5 += 439. / 216. * du_1 * dt - 8. * du_2 * dt + 3680. / 513. * du_3 * dt - 845. / 4104. * du_4 * dt;
                Eigen::VectorXd dx_5;
                Eigen::VectorXd du_5;
                gradient_x(x_5, u_5, dx_5, rho);
                gradient_u(x_5, du_5, alpha);

                Eigen::VectorXd x_6 = this->solver.sol.prime.variables.now;
                x_6 += -8. / 27. * dx_1 * dt + 2. * dx_2 * dt - 3544. / 2565. * dx_3 * dt + 1859. / 4104. * dx_4 * dt - 11. / 40. * dx_5 * dt;
                var_constrained(x_6);
                Eigen::VectorXd u_6 = this->solver.sol.dual.variables.now;
                u_6 += -8. / 27. * du_1 * dt + 2. * du_2 * dt - 3544. / 2565. * du_3 * dt + 1859. / 4104. * du_4 * dt - 11. / 40. * du_5 * dt;
                Eigen::VectorXd dx_6;
                gradient_x(x_6, u_6, dx_6, rho);

                Eigen::VectorXd dh = dt * (16. / 135. * dx_1 + 6656. / 12825. * dx_3 + 28561. / 56430. * dx_4 - 9. / 50. * dx_5 + 2. / 55. * dx_6);
//                Eigen::VectorXd dl = dt * (25. / 216. * dx_1 + 1408. / 2565. * dx_3 + 2197. / 4104. * dx_4 - 1. / 5. * dx_5);
//                double truncate_error = (dl - dh).array().abs().maxCoeff();
//                if(truncate_error > pow(tol_dual, 2.)){
//                    dt *= .9 * pow(pow(tol_dual, 2.) / truncate_error, 1. / 5.);
//                }
//                else if(truncate_error < 1E-4 * pow(tol_dual, 2.)){
//                    dt *= 1.05;
//                }

                this->solver.sol.prime.variables.now += dh;
                var_constrained(this->solver.sol.prime.variables.now);
                Eigen::VectorXd du;
                gradient_u(this->solver.sol.prime.variables.now, du, alpha);
                this->solver.sol.dual.variables.now += du * dt;

                if(loop > 0 && loop % 5 == 0){
                    // Root solving algorithm for exact price
                    solve_iteration_subproblem(rho);

                    // Check convergence for whole problem
                    this->solver.sol.prime.error = this->solver.constraint.Matrix_main * this->solver.sol.prime.variables.now - this->solver.constraint.boundary;
                    this->solver.sol.dual.error = this->solver.sol.prime.price_margin;
                    this->solver.sol.dual.error += rho * this->solver.constraint.Matrix_main.transpose() * this->solver.sol.dual.variables.now;
                    double prime_error_norm = this->solver.sol.prime.error.norm();
                    prime_error_norm /= this->solver.sol.prime.error.size();
                    double dual_error_norm = this->solver.sol.dual.error.norm();
                    dual_error_norm /= this->solver.sol.dual.error.size();
                    if(prime_error_norm < tol_prime && dual_error_norm < tol_dual){
                        break;
                    }

                    if(print_flag){
                        int sys_show = (int) 7. - log(this->statistic.num_variable) / log(10.);
                        sys_show = std::max(sys_show, 0);
                        if(loop % (int) pow(10., sys_show) == 0){
                            std::cout << "Loop:\t" << loop << "\n";
                            std::cout << "Prime Error:\t" << prime_error_norm << "\n";
                            std::cout << "Dual Error:\t" << dual_error_norm << "\n";
                            std::cout << "\n\n";
                        }
                    }

//                    break;
                }

                this->solver.sol.prime.variables.prev = this->solver.sol.prime.variables.now;
                this->solver.sol.dual.variables.prev = this->solver.sol.dual.variables.now;
                loop += 1;
            }

            if(print_flag){
                std::cout << "Total loop:\t" << loop << "\n";
                std::cout << "Prime Error:\t" << this->solver.sol.prime.error.norm() / this->solver.sol.prime.error.size() << "\n";
                std::cout << "Dual Error:\t" << this->solver.sol.dual.error.norm() / this->solver.sol.dual.error.size() << "\n";
                std::cout << "Solution:\n";
                std::cout << (this->obj.transformation.scale.array() * this->solver.sol.prime.variables.now.array() + this->obj.transformation.shift.array()).segment(this->statistic.num_node, this->statistic.num_node).transpose() << "\n";
            }
        }
    };

    // Functions
    void radial_line_problem_set(opf_struct&, int, int, std::complex<double>, double, double, double);
    double eigen_value(Eigen::SparseMatrix <double>);
}