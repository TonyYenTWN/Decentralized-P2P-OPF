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
                void moc_set(){
                    int num_row = 2 * (this->demand.price.size() + this->supply.price.size() - 3);
                    auto price_vec = Eigen::VectorXd(num_row);
                    auto quantity_vec = Eigen::VectorXd(num_row);
                    auto obj_vec = Eigen::VectorXd(num_row);

                    // Initialization for the loop
                    int price_demand_ID = 0;
                    int price_supply_ID = 0;
                    int price_ID = 0;
                    double obj_value = -(this->demand.price.segment(1, this->demand.price.size() - 2).array() * this->demand.quantity.segment(1, this->demand.quantity.size() - 2).array()).sum();

                    // First item when price = price_min
                    price_vec(0) = -std::numeric_limits<double>::infinity();
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
                    price_vec(price_ID) = std::numeric_limits<double>::infinity();
                    quantity_vec(price_ID) = quantity_vec(price_ID - 1);
                    obj_vec(price_ID) = obj_value;

                    // Extend the domain to (-inf, inf)
                    this->moc.price = price_vec;
                    this->moc.quantity = quantity_vec;
                    this->moc.obj = obj_vec;
                }
            };

            std::vector <cost_func_struct> cost_funcs;
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

                opf_struct::obj_struct::cost_func_struct cost_func;
                cost_func.moc.price = Eigen::VectorXd(4);
                cost_func.moc.quantity = Eigen::VectorXd(4);
                cost_func.moc.obj = Eigen::VectorXd::Zero(4);

                cost_func.moc.price << -std::numeric_limits<double>::infinity(), 0., 0., std::numeric_limits<double>::infinity();
                cost_func.moc.quantity << -current_limit, -current_limit, current_limit, current_limit;

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

        // Calculate objective value
        void get_obj(){
            this->solver.sol.obj_value = 0.;
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                // Bisection method for finding location of solution
                int quan_min_ID = 0;
                int quan_max_ID = this->obj.cost_funcs[var_iter].moc.price.size() - 1;
                Eigen::Vector2i gap_price_ID(quan_min_ID, quan_max_ID);
                int mid_price_ID = gap_price_ID.sum() / 2;
                Eigen::Vector2d gap_rent;
                while(gap_price_ID(1) - gap_price_ID(0) > 1){
                    gap_rent(0) = this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(0)) - this->solver.sol.variable.value.now(var_iter);
                    gap_rent(1) = this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(1)) - this->solver.sol.variable.value.now(var_iter);

                    double mid_rent;
                    mid_rent = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID) - this->solver.sol.variable.value.now(var_iter);

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
        }

        // Subsolver for price
        // DP: min g*(z) - <z - z_{t - 1}, A %*% (2 * x_t - x_{t - 1})> + .5 * ||z - z_{t - 1}||^2_{M_2}
        // KKT: q(z) - A %*% (2 * x_t - x_{t - 1}) + 1 / rho * A %*% A^T %*% (z - z_{t - 1}) = 0
        void solve_iteration_subproblem(double rho){
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                double inner_prod = 0.;

                // add from PSD
                int entry_iter = 0;
                while(entry_iter < this->solver.constraint.PSD_terms.row[var_iter].size()){
                    auto pair_temp = this->solver.constraint.PSD_terms.row[var_iter][entry_iter];
                    if(pair_temp.first != var_iter){
                        inner_prod += pair_temp.second * this->solver.sol.price.value.now(pair_temp.first);
                    }
                    inner_prod -= pair_temp.second * this->solver.sol.price.value.prev(pair_temp.first);
                    entry_iter += 1;
                }

                inner_prod /= rho;
                inner_prod -= 2 * this->solver.sol.variable.value.now(var_iter) - this->solver.sol.variable.value.prev(var_iter);

                // Bisection method for finding KKT point
                int quan_min_ID = 0;
                int quan_max_ID = this->obj.cost_funcs[var_iter].moc.price.size() - 1;
                Eigen::Vector2i gap_price_ID(quan_min_ID, quan_max_ID);
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

                // Check binding term
                if(gap_price_ID(0) % 2 == 1){
                    // price fixed case
                    this->solver.sol.price.value.now(var_iter) = this->obj.cost_funcs[var_iter].moc.price(mid_price_ID);
                    this->solver.sol.price.grad_margin(var_iter) = -this->solver.sol.price.value.now(var_iter);
                    this->solver.sol.price.grad_margin(var_iter) *= this->solver.constraint.sensitivity_var(var_iter) / rho;
                    this->solver.sol.price.grad_margin(var_iter) -= inner_prod;
                }
                else{
                    // quantity fixed case
                    this->solver.sol.price.grad_margin(var_iter) = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID);
                    this->solver.sol.price.value.now(var_iter) = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID);
                    this->solver.sol.price.value.now(var_iter) += inner_prod;
                    this->solver.sol.price.value.now(var_iter) *= -rho / this->solver.constraint.sensitivity_var(var_iter);
                }
            }
        }

        // Main solver
        void solve_root(double tol_prime, double tol_dual, double theta_limit, bool print_flag = 1){
            // Initialization
            double rho = 100.;
            this->solver.sol.state.value.prev = Eigen::VectorXd::Zero(this->statistic.num_state);
            this->solver.sol.state.value.now = Eigen::VectorXd::Zero(this->statistic.num_state);
            this->solver.sol.variable.value.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.variable.value.now = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.price.value.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.price.value.now = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.price.grad_margin = Eigen::VectorXd::Zero(this->statistic.num_variable);

            // Main loop
            // PP: min f(x) + <x - x_{t - 1}, A^T * z_{t - 1}> + .5 * || x - x_{t - 1} ||^2_{M_1}
            // KKT: p(x) + A^T %*% z_{t - 1} + rho * (x - x_{t - 1}) = 0
            // DP: min g*(z) - <z - z_{t - 1}, A %*% (2 * x_t - x_{t - 1})> + .5 * ||z - z_{t - 1}||^2_{M_2}
            // KKT: q(z) - A %*% (2 * x_t - x_{t - 1}) + 1 / rho * A %*% A^T %*% (z - z_{t - 1}) = 0
            int loop = 0;
            while(true){
                // Update states
                this->solver.sol.state.value.now -= 1. / rho * this->solver.constraint.Mat_main.transpose() * this->solver.sol.price.value.now;
                this->solver.sol.state.value.now = this->solver.sol.state.value.now.array().min(theta_limit);
                this->solver.sol.state.value.now = this->solver.sol.state.value.now.array().max(-theta_limit);

                // Update variables
                this->solver.sol.variable.value.now = this->solver.constraint.Mat_main * this->solver.sol.state.value.now;

                // Update prices
                int sub_loop = 0;
                while(true){
                    Eigen::VectorXd price_prev_temp = this->solver.sol.price.value.now;
                    solve_iteration_subproblem(rho);

                    // Check convergence for subproblem
                    double tol_sub = 1E-2 * tol_dual;
                    Eigen::VectorXd dual_error = this->solver.sol.price.value.now - price_prev_temp;
                    dual_error = this->solver.constraint.PSD_main * dual_error;
                    dual_error /= rho;
                    if(dual_error.norm() / this->solver.sol.price.value.now.size() < tol_sub){
                        break;
                    }

                    if(sub_loop > -1){
                        break;
                    }
                    sub_loop += 1;
                }

                // Check convergence for whole problem
                this->solver.sol.state.error = this->solver.sol.state.value.now - this->solver.sol.state.value.prev;
                this->solver.sol.price.error = this->solver.constraint.PSD_main * (this->solver.sol.price.value.now - this->solver.sol.price.value.prev);
                double prime_error_norm = this->solver.sol.state.error.norm() * rho;
                prime_error_norm /= this->solver.sol.state.error.size();
                double dual_error_norm = this->solver.sol.price.error.norm() / rho;
                dual_error_norm /= this->solver.sol.price.error.size();
                if(prime_error_norm < tol_prime && dual_error_norm < tol_dual){
                    break;
                }

                // Update previous values
                this->solver.sol.state.value.prev = this->solver.sol.state.value.now;
                this->solver.sol.variable.value.prev = this->solver.sol.variable.value.now;
                this->solver.sol.price.value.prev = this->solver.sol.price.value.now;

//                // Update rho
//                double drho = this->solver.sol.state.error.array().abs().maxCoeff();
//                drho /= this->solver.sol.price.error.array().abs().maxCoeff();
//                drho = rho * pow(drh)
//                rho *= pow(drho, .5);
//                rho = std::max(rho, 1.);
//                rho = std::min(rho, 1E6);

                // Print progress
                if(print_flag){
                    int sys_show = (int) 7. - log(this->statistic.num_variable) / log(10.);
                    sys_show = std::max(sys_show, 0);
                    if(loop % (int) pow(10., sys_show) == 0){
                        std::cout << "Loop:\t" << loop << "\n";
                        std::cout << "Rho:\t" << rho << "\n";
                        std::cout << "Prime Error:\t" << prime_error_norm << "\n";
                        std::cout << "Dual Error:\t" << dual_error_norm << "\n";
                        std::cout << "\n\n";
                    }
                }
                loop += 1;
            }

            // Calculate objective value
            get_obj();

            if(print_flag){
                std::cout << "Total loop:\t" << loop << "\n";
                std::cout << "Prime Error:\t" << this->solver.sol.state.error.norm() / this->solver.sol.state.error.size() << "\n";
                std::cout << "Dual Error:\t" << this->solver.sol.price.error.norm() / this->solver.sol.price.error.size() << "\n";
                std::cout << "Objective value:\t" << this->solver.sol.obj_value << "\n";
                std::cout << "Solution:\n";
                std::cout << (this->solver.sol.variable.value.now.head(this->statistic.num_node)).transpose() << "\n";
                std::cout << "\n";
           }
        }
    };

    // Functions
    void radial_line_problem_set(opf_struct&, int, int, std::complex<double>, double, double, double, double);
}
