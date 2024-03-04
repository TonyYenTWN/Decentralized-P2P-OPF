// header for ADMM

#pragma once

// STL
#include <complex>
#include <iostream>
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

    struct SparseMatrix{
        std::vector <std::vector <std::pair <int, double>>> row;
        std::vector <std::vector <std::pair <int, double>>> col;

        void dim_set(int num_row, int num_col){
            this->row = std::vector <std::vector <std::pair <int, double>>> (num_row);
            this->col = std::vector <std::vector <std::pair <int, double>>> (num_col);
        }

        void read_from_Dense(Eigen::MatrixXd Dense){
            // Initialization
            for(int row_iter = 0; row_iter < Dense.rows(); ++ row_iter){
                this->row[row_iter].reserve(Dense.cols());
            }
            for(int col_iter = 0; col_iter < Dense.cols(); ++ col_iter){
                this->col[col_iter].reserve(Dense.rows());
            }

            // Scan through the dense matrix
            for(int row_iter = 0; row_iter < Dense.rows(); ++ row_iter){
                for(int col_iter = 0; col_iter < Dense.cols(); ++ col_iter){
                    if(Dense(row_iter, col_iter) != 0.){
                        double coeff = Dense(row_iter, col_iter);
                        this->row[row_iter].push_back(std::pair <int, double> (col_iter, coeff));
                        this->col[col_iter].push_back(std::pair <int, double> (row_iter, coeff));
                    }
                }
            }
        }
    };

    struct opf_struct{
        // System statistic
        struct statistic_struct{
            int num_node;
            int num_line;
            int num_state;
            int num_variable;
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

                // Function to set merit order curve
                void moc_set(Eigen::Vector2d price_range){
                    int num_row = 2 * (this->demand.price.size() + this->supply.price.size() - 3);
                    auto price_vec = Eigen::VectorXd(num_row);
                    auto quantity_vec = Eigen::VectorXd(num_row);
                    auto obj_vec = Eigen::VectorXd(num_row);
                    this->moc.price = Eigen::VectorXd(num_row + 2);
                    this->moc.quantity = Eigen::VectorXd(num_row + 2);
                    this->moc.obj = Eigen::VectorXd(num_row + 2);

                    // Initialization for the loop
                    int price_demand_ID = 0;
                    int price_supply_ID = 0;
                    int price_ID = 0;
                    double obj_value = -(this->demand.price.segment(1, this->demand.price.size() - 2).array() * this->demand.quantity.segment(1, this->demand.quantity.size() - 2).array()).sum();

                    // First item when price = price_min
                    price_vec(0) = price_range(0);
                    quantity_vec(0) = -this->demand.quantity.sum();
                    obj_vec(0) = obj_value;

                    while(price_ID < num_row - 2){
                        price_ID += 1;

                        // Next demand marginal price lower than supply; add that demand quantity to residual load
                        if(this->demand.price(price_demand_ID + 1) <= this->supply.price(price_supply_ID + 1)){
                            price_demand_ID += price_demand_ID < this->demand.price.size() - 2;
                            price_vec(price_ID) = this->demand.price(price_demand_ID);
                            quantity_vec(price_ID) = quantity_vec(price_ID - 1);
                            obj_vec(price_ID, 2) = obj_value;

                            price_ID += 1;
                            obj_value += this->demand.price(price_demand_ID) * this->demand.quantity(price_demand_ID);
                            price_vec(price_ID) = this->demand.price(price_demand_ID);
                            quantity_vec(price_ID) = quantity_vec(price_ID - 1) + this->demand.quantity(price_demand_ID);
                            obj_vec(price_ID) = obj_value;
                        }
                        // Next supply marginal price lower than demand; add that supply quantity to residual load
                        else{
                            price_supply_ID += price_supply_ID < this->supply.price.size() - 2;
                            price_vec(price_ID) = this->supply.price(price_supply_ID);
                            quantity_vec(price_ID) = quantity_vec(price_ID - 1);
                            obj_vec(price_ID, 2) = obj_value;

                            price_ID += 1;
                            obj_value += this->supply.price(price_supply_ID) * this->supply.quantity(price_supply_ID);
                            price_vec(price_ID) = this->supply.price(price_supply_ID);
                            quantity_vec(price_ID) = quantity_vec(price_ID - 1) + this->supply.quantity(price_supply_ID);
                            obj_vec(price_ID) = obj_value;
                        }
                    }

                    // Last item when price = max
                    price_ID += 1;
                    price_vec(price_ID) = price_range(1);
                    quantity_vec(price_ID) = quantity_vec(price_ID - 1);
                    obj_vec(price_ID) = obj_value;

                    // Extend the domain to (-inf, inf)
                    this->moc.price << price_range(0), price_vec, price_range(1);
                    this->moc.quantity << -std::numeric_limits<double>::infinity(), quantity_vec, std::numeric_limits<double>::infinity();
                    this->moc.obj << std::numeric_limits<double>::infinity(), obj_vec, std::numeric_limits<double>::infinity();
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
            // Augmented Lagrangian with Prime-Dual hybrid gradient:
            // min f(x) + g(A %*% x)
            // PD: min f(x) + <x - x_{t - 1}, A^T * z_{t - 1}> + .5 * || x - x_{t - 1} ||^2_{M_1}
            // DD: min g*(z) - <z - z_{t - 1}, A %*% (2 * x_t - x_{t - 1})> + .5 * ||z - z_{t - 1}||^2_{M_2}
            // x: voltages at each node
            // A: equality constraints transforming voltages to nodal power sources and line currents
            // g*(z): convex conjugate of g(x)
            // g*(z) = sup{z^T * A %*% x - g(A * x)}
            // (g*(z))' = x(z) (See notes)
            // M_1: pre-conditioner for prime problem
            // M_1 = rho * I, rho > 1
            // M_2: pre-conditioner for dual problem
            // M_2 = 1 / rho * A %*% A^T

            struct constraint_struct{
                // Main matrix
                Eigen::SparseMatrix <double> Mat_main;
                Eigen::SparseMatrix <double> PSD_main;
                SparseMatrix Mat_terms;
                SparseMatrix PSD_terms;
                Eigen::VectorXd sensitivity_var;
                Eigen::VectorXd boundary;
            };
            constraint_struct constraint;

            // Solution
            struct sol_struct{
                struct variable_struct{
                    struct value_struct{
                        Eigen::VectorXd now;
                        Eigen::VectorXd prev;
                    };
                    value_struct value;

                    Eigen::VectorXd grad_margin;
                    Eigen::VectorXd error;
                };

                variable_struct state;
                variable_struct variable;
                variable_struct price;
                double obj_value;
            };
            sol_struct sol;
        };
        solver_struct solver;

        // Functions
        // Set size for the LP
        void problem_size_parmeters_set(int factor = 1){
            this->statistic.num_state = this->statistic.num_node;
            this->statistic.num_variable = factor * (this->statistic.num_node + this->statistic.num_line);
        }

        // Initialize merit order curves and set values for voltage and current
        void moc_initialize(double current_limit){
            this->obj.cost_funcs = std::vector <obj_struct::cost_func_struct> (this->statistic.num_variable);

            // Line current boundaries
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = this->statistic.num_node + line_iter;
                double price_range = this->obj.price_range(1) - this->obj.price_range(0);

                opf_struct::obj_struct::cost_func_struct cost_func;
                cost_func.moc.price = Eigen::VectorXd(6);
                cost_func.moc.quantity = Eigen::VectorXd(6);
                cost_func.moc.obj = Eigen::VectorXd::Zero(6);

                cost_func.moc.price << -price_range, -price_range, 0., 0., price_range, price_range;
                cost_func.moc.quantity << -std::numeric_limits<double>::infinity(), -current_limit, -current_limit, current_limit, current_limit, std::numeric_limits<double>::infinity();
                cost_func.moc.obj(0) = std::numeric_limits<double>::infinity();
                cost_func.moc.obj(5) = std::numeric_limits<double>::infinity();

                this->obj.cost_funcs[var_ID] = cost_func;
            }
        }

        // Main matrix initialization for DC OPF
        void DC_Matrix_main_set(){
            // States: V
            // Variables: S, I
            // Equalities: Nodal power flow, line current
            this->solver.constraint.Mat_terms.dim_set(this->statistic.num_variable, this->statistic.num_state);
            this->solver.constraint.PSD_terms.dim_set(this->statistic.num_variable, this->statistic.num_variable);
            this->solver.constraint.sensitivity_var = Eigen::VectorXd (this->statistic.num_variable);

            // Set the main matrix
            Eigen::MatrixXd Mat = Eigen::MatrixXd::Zero(this->statistic.num_variable, this->statistic.num_state);

            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int from_ID = this->network.topology[line_iter](0);
                int to_ID = this->network.topology[line_iter](1);
                double y_l = this->network.line_conductance(line_iter).imag();

                // Nodal power flow
                // Y %*% V = S
                // Non-diagonal terms
                Mat(from_ID, to_ID) -= y_l;
                Mat(to_ID, from_ID) -= y_l;

                // Diagonal terms
                Mat(to_ID, to_ID) += y_l;
                Mat(from_ID, from_ID) += y_l;

                // Line Current Equation
                // Y_l %*% NL %*% V = I
                int I_ID = this->statistic.num_node + line_iter;
                Mat(I_ID, from_ID) = y_l;
                Mat(I_ID, to_ID) = -y_l;
            }
            Eigen::MatrixXd PSD = Mat * Mat.transpose();

            // Record the diagonal terms of PSD
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                this->solver.constraint.sensitivity_var(var_iter) = PSD(var_iter, var_iter);
            }

            // Turn dense to sparse matrix
            this->solver.constraint.Mat_main = Mat.sparseView();
            this->solver.constraint.PSD_main = PSD.sparseView();

            // Record non-zero terms for matrices
            this->solver.constraint.Mat_terms.read_from_Dense(Mat);
            this->solver.constraint.PSD_terms.read_from_Dense(PSD);
        }

        void solve_root(double tol_prime, double tol_dual, double theta_limit, bool print_flag = 1){
            // Initialization
            double rho = 1.;
            this->solver.sol.state.value.prev = Eigen::VectorXd::Zero(this->statistic.num_state);
            this->solver.sol.state.value.now = Eigen::VectorXd::Zero(this->statistic.num_state);
            this->solver.sol.variable.value.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.variable.value.now = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.price.value.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.price.value.now = Eigen::VectorXd::Zero(this->statistic.num_variable);

            // Main loop
            // PP: min f(x) + <x - x_{t - 1}, A^T * z_{t - 1}> + .5 * || x - x_{t - 1} ||^2_{M_1}
            // KKT: p(x) + A^T %*% z_{t - 1} + rho * (x - x_{t - 1}) = 0
            // DP: min g*(z) - <z - z_{t - 1}, A %*% (2 * x_t - x_{t - 1})> + .5 * ||z - z_{t - 1}||^2_{M_2}
            // KKT: q(z) - A %*% (2 * x_t - x_{t - 1}) + 1 / rho * A %*% A^T %*% (z - z_{t - 1}) = 0
            int loop = 0;
            while(true){
                // Update prices
                for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                    double inner_prod = 0.;

                    // add from PSD
                    int entry_iter = 0;
                    while(entry_iter < this->solver.constraint.PSD_terms.row[var_iter].size()){
                        auto pair_temp = this->solver.constraint.PSD_terms.row[var_iter][entry_iter];
                        if(pair_temp.first != var_iter){
                            inner_prod += pair_temp.second * this->solver.sol.price.value.now(entry_iter);
                        }
                        inner_prod -= pair_temp.second * this->solver.sol.price.value.prev(entry_iter);
                        entry_iter += 1;
                    }

                    inner_prod /= rho;
                    inner_prod -= 2 * this->solver.sol.variable.value.now(var_iter) - this->solver.sol.variable.value.prev(var_iter);

                    // Bisection method for finding KKT point
                    Eigen::Vector2i gap_price_ID(0, this->obj.cost_funcs[var_iter].moc.price.size() - 1);
                    int mid_price_ID = gap_price_ID.sum() / 2;
                    while(gap_price_ID(1) - gap_price_ID(0) > 1){
                        Eigen::Vector2d gap_rent;
                        gap_rent(0) = this->obj.cost_funcs[var_iter].moc.price(gap_price_ID(0));
                        gap_rent(0) *= this->solver.constraint.sensitivity_var(var_iter) / rho;
                        gap_rent(0) += this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(0));
                        gap_rent(0) += inner_prod;
                        gap_rent(1) = this->obj.cost_funcs[var_iter].moc.price(gap_price_ID(1));
                        gap_rent(1) *= this->solver.constraint.sensitivity_var(var_iter) / rho;
                        gap_rent(1) += this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(1));
                        gap_rent(1) += inner_prod;

                        double mid_rent;
                        mid_rent = this->obj.cost_funcs[var_iter].moc.price(mid_price_ID);
                        mid_rent *= this->solver.constraint.sensitivity_var(var_iter) / rho;
                        mid_rent += this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID);
                        mid_rent += inner_prod;

                        // Bisection root-search method
                        if(mid_rent * gap_rent(0) <= 0){
                            gap_price_ID(1) = mid_price_ID;
                        }
                        else{
                            gap_price_ID(0) = mid_price_ID;
                        }

                        mid_price_ID = gap_price_ID.sum() / 2;
                    }
                }

                // Update states
                this->solver.sol.state.value.now -= 1. / rho * this->solver.constraint.Mat_main.transpose() * this->solver.sol.price.value.now;
                this->solver.sol.state.value.now = this->solver.sol.state.value.now.array().min(theta_limit);
                this->solver.sol.state.value.now = this->solver.sol.state.value.now.array().max(-theta_limit);

                // Update variables
                this->solver.sol.variable.value.now = this->solver.constraint.Mat_main * this->solver.sol.state.value.now;

                break;
            }
        }
    };

    // Functions
    void radial_line_problem_set(opf_struct&, int, int, std::complex<double>, double, double, double, double);
}
