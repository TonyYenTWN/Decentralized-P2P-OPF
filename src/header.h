// header for ADMM

#pragma once

// STL
#include <complex>
#include <iostream>
#include <map>
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
            // Subset information (for lower level networks)
            struct network_map_struct{
                int current_level;
                std::map <int, int> connection_nodes; // argument: current level node ID; mapped value: upper level connected node ID (-1 if none)
            };
            network_map_struct network_map;

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
                    Eigen::VectorXd gap;
                    Eigen::VectorXd obj;
                };
                moc_struct moc;

                // Function to merit order curve
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

            Eigen::Vector2d price_range;
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
                std::vector <std::vector <std::pair <int, double>>> Mat_main_terms;
                std::vector <std::vector <std::pair <int, double>>> PSD_main_terms;
                Eigen::VectorXd sensitivity_var;
                Eigen::VectorXd boundary_0;
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

        // Linear transformation of the variables and objective coefficients
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
        }

        // Initialize merit order curves and set values for voltage and current
        void moc_initialize(double theta_limit, double current_limit, double penalty_price_voltage){
            this->obj.cost_funcs = std::vector <obj_struct::cost_func_struct> (this->statistic.num_variable);
            this->obj.price_range << -500., 3000.;
            double price_gap = 1.05 * (this->obj.price_range(1) - this->obj.price_range(0));
            double quantity_inflex = 1000.;

            // Line current boundaries
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = 2 * this->statistic.num_node + line_iter;

                // Set bid functions for suuply
                this->obj.cost_funcs[var_ID].supply.price = Eigen::VectorXd(4);
                this->obj.cost_funcs[var_ID].supply.quantity = Eigen::VectorXd::Zero(4);
                this->obj.cost_funcs[var_ID].supply.price << -std::numeric_limits<double>::infinity(), 0., price_gap, std::numeric_limits<double>::infinity();
                this->obj.cost_funcs[var_ID].supply.quantity << 0., current_limit, quantity_inflex, 0.;

                // Set bid functions for demand
                this->obj.cost_funcs[var_ID].demand.price = Eigen::VectorXd(4);
                this->obj.cost_funcs[var_ID].demand.quantity = Eigen::VectorXd(4);
                this->obj.cost_funcs[var_ID].demand.price << -std::numeric_limits<double>::infinity(), -price_gap, 0., std::numeric_limits<double>::infinity();
                this->obj.cost_funcs[var_ID].demand.quantity << 0., quantity_inflex, current_limit, 0.;

                // Set merit order curve for residual load
                this->obj.cost_funcs[var_ID].moc_set();
            }

            if(penalty_price_voltage == std::numeric_limits<double>::infinity()){
                // Phase angle boundaries
                for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                    int var_ID = node_iter;
                    opf_struct::obj_struct::cost_func_struct cost_func;
                    cost_func.moc.price = Eigen::VectorXd(4);
                    cost_func.moc.quantity = Eigen::VectorXd(4);
                    cost_func.moc.obj = Eigen::VectorXd::Zero(4);

                    cost_func.moc.price << -std::numeric_limits<double>::infinity(), 0., 0., std::numeric_limits<double>::infinity();
                    cost_func.moc.quantity << -theta_limit, -theta_limit, theta_limit, theta_limit;

                    this->obj.cost_funcs[var_ID] = cost_func;
                }

                return;
            }

            // Phase angle boundaries
            for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                int var_ID = node_iter;

                // Set bid functions for suuply
                this->obj.cost_funcs[var_ID].supply.price = Eigen::VectorXd(4);
                this->obj.cost_funcs[var_ID].supply.quantity = Eigen::VectorXd::Zero(4);
                this->obj.cost_funcs[var_ID].supply.price << -std::numeric_limits<double>::infinity(), 0., penalty_price_voltage, std::numeric_limits<double>::infinity();
                this->obj.cost_funcs[var_ID].supply.quantity << 0., theta_limit, quantity_inflex, 0.;

                // Set bid functions for demand
                this->obj.cost_funcs[var_ID].demand.price = Eigen::VectorXd(4);
                this->obj.cost_funcs[var_ID].demand.quantity = Eigen::VectorXd(4);
                this->obj.cost_funcs[var_ID].demand.price << -std::numeric_limits<double>::infinity(), -penalty_price_voltage, 0., std::numeric_limits<double>::infinity();
                this->obj.cost_funcs[var_ID].demand.quantity << 0., quantity_inflex, theta_limit, 0.;

                // Set merit order curve for residual load
                this->obj.cost_funcs[var_ID].moc_set();
            }
        }

        // Main matrix initialization for DC OPF
        void DC_Matrix_main_set(){
            // Variables: V, S, I
            // Constraints: Node balance, line current
            this->solver.constraint.Mat_main_terms = std::vector <std::vector <std::pair <int, double>>> (this->statistic.num_variable);
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->solver.constraint.Mat_main_terms[var_iter].reserve(this->statistic.num_constraint);
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

            // Set boundary for equality constraints; transformation considered
            // Ax = c
            // A(mx' + x_0) = c
            // [A][m]x' = c - A * x_0
            this->solver.constraint.boundary_0 = -Mat * this->obj.transformation.shift;
            this->solver.constraint.boundary = this->solver.constraint.boundary_0;

            // Apply transformation to the matrix
            Eigen::MatrixXd Diag = Eigen::MatrixXd::Zero(this->statistic.num_variable, this->statistic.num_variable);
            Eigen::MatrixXd Diag_inv = Eigen::MatrixXd::Zero(this->statistic.num_variable, this->statistic.num_variable);
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                Diag(var_iter, var_iter) = this->obj.transformation.scale(var_iter);
                Diag_inv(var_iter, var_iter) = 1. / Diag(var_iter, var_iter);
            }
            Mat = Mat * Diag;
            Eigen::MatrixXd PSD = Mat.transpose() * Mat;

            // Turn dense to sparse matrix
            this->solver.constraint.Matrix_main = Mat.sparseView();
            this->solver.constraint.PSD_main = PSD.sparseView();

            // Record non-zero terms for Mat
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                for(int constr_iter = 0; constr_iter < this->statistic.num_constraint; ++ constr_iter){
                    if(Mat(constr_iter, var_iter) != 0.){
                        this->solver.constraint.Mat_main_terms[var_iter].push_back(std::pair <int, double> (constr_iter, Mat(constr_iter, var_iter)));
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
        }

        // Update boundary terms
        // Only works for bisplit now!!
        void boundary_update(double boundary_current = 0., bool right_end = 1){
            int node_ID = right_end * (this->statistic.num_node - 1);
            this->solver.constraint.boundary(node_ID) = this->solver.constraint.boundary_0(node_ID);
            this->solver.constraint.boundary(node_ID) += boundary_current * (2 * right_end - 1);
        }

        // Price gap set (rho must be fixed!!)
        void price_gap_set(double rho){
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->obj.cost_funcs[var_iter].moc.gap = this->obj.cost_funcs[var_iter].moc.price / rho;
                this->obj.cost_funcs[var_iter].moc.gap += this->solver.constraint.sensitivity_var(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity;
            }
        }

        // Solver function: direct search the root
        void solve_iteration_subproblem(double rho, Eigen::VectorXd &x, Eigen::VectorXd &u){
            // Update primary variables and price margin
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                double inner_prod = 0.;
                // add from PSD
                int entry_iter = 0;
                while(entry_iter < this->solver.constraint.PSD_main_terms[var_iter].size()){
                    int row_ID = this->solver.constraint.PSD_main_terms[var_iter][entry_iter].first;
                    double coeff = this->solver.constraint.PSD_main_terms[var_iter][entry_iter].second;
                    double value = coeff * x(row_ID);
                    inner_prod += value;
                    entry_iter += 1;
                }

                // add from Mat
                entry_iter = 0;
                while(entry_iter < this->solver.constraint.Mat_main_terms[var_iter].size()){
                    int constr_ID = this->solver.constraint.Mat_main_terms[var_iter][entry_iter].first;
                    double coeff = this->solver.constraint.Mat_main_terms[var_iter][entry_iter].second;
                    double value = coeff * (u(constr_ID) - this->solver.constraint.boundary(constr_ID));
                    inner_prod += value;
                    entry_iter += 1;
                }

                // Bisection method for finding KKT point
                Eigen::Vector2i gap_price_ID(0, this->obj.cost_funcs[var_iter].moc.price.size() - 1);
                int mid_price_ID = gap_price_ID.sum() / 2;
                while(gap_price_ID(1) - gap_price_ID(0) > 1){
                    Eigen::Vector2d gap_rent;
                    gap_rent(0) = this->obj.cost_funcs[var_iter].moc.gap(gap_price_ID(0)) + inner_prod;
                    gap_rent(1) = this->obj.cost_funcs[var_iter].moc.gap(gap_price_ID(1)) + inner_prod;

                    double mid_rent;
                    mid_rent = this->obj.cost_funcs[var_iter].moc.gap(mid_price_ID) + inner_prod;

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
                    x(var_iter) = -this->solver.sol.prime.price_margin(var_iter) / rho;
                    x(var_iter) -= inner_prod;
                    x(var_iter) /= this->solver.constraint.sensitivity_var(var_iter);
                }
                else{
                    // quantity fixed case
                    x(var_iter) = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID);
                    this->solver.sol.prime.price_margin(var_iter) = -rho * (this->solver.constraint.sensitivity_var(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID) + inner_prod);
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
            price_gap_set(rho);

            // Main loop
            // f(x) + .5 * rho * ||A %*% x - c + u||^2
            // f(x) + .5 * rho * (M_{ij} x_i x_j + 2 * A_{ij} b_i x_j + ...)
            // KKT: p_i + rho * (M_{ij} x_j + t(A)_{ij} b_j) = 0
            int loop = 0;
            while(true){
                // Update prime variables and associate prices
                int sub_loop = 0;
                while(true){
                    solve_iteration_subproblem(rho, this->solver.sol.prime.variables.now, this->solver.sol.dual.variables.now);

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
                    if(loop % (int) pow(10., sys_show) == 0){
                        std::cout << "Loop:\t" << loop << "\n";
                        std::cout << "Prime Error:\t" << prime_error_norm << "\n";
                        std::cout << "Dual Error:\t" << dual_error_norm << "\n";
                        std::cout << "\n\n";
                    }
                }
                loop += 1;
            }

            // Calculate objective value
            this->solver.sol.obj_value = 0.;
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                // Bisection method for finding location of solution
                Eigen::Vector2i gap_price_ID(0, this->obj.cost_funcs[var_iter].moc.price.size() - 1);
                int mid_price_ID = gap_price_ID.sum() / 2;
                Eigen::Vector2d gap_rent;
                while(gap_price_ID(1) - gap_price_ID(0) > 1){
                    gap_rent(0) = this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(0)) - this->solver.sol.prime.variables.now(var_iter);
                    gap_rent(1) = this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(1)) - this->solver.sol.prime.variables.now(var_iter);

                    double mid_rent;
                    mid_rent = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID) - this->solver.sol.prime.variables.now(var_iter);

                    // Bisection root-search method
                    if(mid_rent * gap_rent(0) <= 0){
                        gap_price_ID(1) = mid_price_ID;
                    }
                    else{
                        gap_price_ID(0) = mid_price_ID;
                    }

                    mid_price_ID = gap_price_ID.sum() / 2;
                }

                double obj_temp = 0.;
                obj_temp += this->obj.cost_funcs[var_iter].moc.obj(gap_price_ID(0)) * gap_rent(1);
                obj_temp -= this->obj.cost_funcs[var_iter].moc.obj(gap_price_ID(1)) * gap_rent(0);
                obj_temp /= gap_rent(1) - gap_rent(0);

                this->solver.sol.obj_value += obj_temp;
            }

            if(print_flag){
                std::cout << "Total loop:\t" << loop << "\n";
                std::cout << "Prime Error:\t" << this->solver.sol.prime.error.norm() / this->solver.sol.prime.error.size() << "\n";
                std::cout << "Dual Error:\t" << this->solver.sol.dual.error.norm() / this->solver.sol.dual.error.size() << "\n";
                std::cout << "Objective value:\t" << this->solver.sol.obj_value << "\n";
                std::cout << "Solution:\n";
                std::cout << (this->obj.transformation.scale.array() * this->solver.sol.prime.variables.now.array() + this->obj.transformation.shift.array()).segment(this->statistic.num_node, this->statistic.num_node).transpose() << "\n";
                std::cout << "\n";
           }
        }

        // Might be used later in parallel
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
    };

    // Functions
    void radial_line_problem_set(opf_struct&, int, int, std::complex<double>, double, double, double, double);
}
