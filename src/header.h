// header for ADMM

#pragma once
// STL
#include <complex>
#include <vector>
#include <iostream>

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

            double theta_mul;
            std::vector <cost_func_struct> cost_funcs;
        };
        obj_struct obj;

        // Solver
        struct solver_struct{
            // Augmented Lagrangian:
            // L(x, u) = f(x) + .5 * \rho * ||A %*% x + u||^2
            // x = scaled voltage (V), power source / sink (S), and line current (I)
            // u = dual variables of the constraints
            // f(x): cost function associated with V, S, and I (f(V) and f(I) are well functions)

            // Main matrix
            Eigen::SparseMatrix <double> Matrix_main;

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

        void theta_mul_calc(){
            double log_y_l = 0;
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                log_y_l += log(std::abs(this->network.line_conductance(line_iter)));
            }
            this->obj.theta_mul = exp(log_y_l / this->statistic.num_line);
        }

        // Main matrix initialization for DC OPF
        void DC_Matrix_main_set(){
            // Variables: V, S, I
            // Constraints: Node balance, line current
            int num_entry = this->statistic.num_node + 5 * this->statistic.num_line;

            // Set the main matrix
            this->solver.Matrix_main = Eigen::SparseMatrix <double> (this->statistic.num_constraint, this->statistic.num_variable);
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
            // set U = mean(Y_l) V
            // 1/mean(Y_l) * Y_l %*% NL %*% U - I = 0
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = 2 * this->statistic.num_node + line_iter;
                int constr_ID = this->statistic.num_node + line_iter;
                double scale = this->network.line_conductance(line_iter).imag() / this->obj.theta_mul;

                Matrix_main_trip.push_back(Eigen::Triplet <double> (constr_ID, this->network.topology[line_iter](0), scale));
                Matrix_main_trip.push_back(Eigen::Triplet <double> (constr_ID, this->network.topology[line_iter](1), -scale));
                Matrix_main_trip.push_back(Eigen::Triplet <double> (constr_ID, var_ID, -1));
            }

            // Put Everything in to the matrix
            this->solver.Matrix_main.setFromTriplets(Matrix_main_trip.begin(), Matrix_main_trip.end());
        }

        // Solver function
        void solve(double tol_prime, double tol_dual){
            // Calculate sensitivity for each variable
            Eigen::VectorXd sensitivity_x(this->statistic.num_variable);
            Eigen::SparseMatrix <double> Diagonal(this->statistic.num_variable, this->statistic.num_variable);
            std::vector <Eigen::Triplet <double>> Diagonal_trip;
            Diagonal_trip.reserve(this->statistic.num_variable);
            for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                Eigen::VectorXd green_vec = Eigen::VectorXd::Zero(this->statistic.num_variable);
                green_vec(var_iter) = 1.;
                Eigen::VectorXd A_col = this->solver.Matrix_main * green_vec;
                sensitivity_x(var_iter) = (A_col.array() * A_col.array()).sum();
                Diagonal_trip.push_back(Eigen::Triplet <double> (var_iter, var_iter, 1. / sensitivity_x(var_iter)));
            }
            Diagonal.setFromTriplets(Diagonal_trip.begin(), Diagonal_trip.end());

            // Calculate appropriate regulation / weight factor for Jacobi method
            // for omega < 2 / (lambda_max + lambda_min) and PSD matrix, Jacobi method converges globally
            // theoretical optimum is omega = 2 / (lambda_max + lambda_min), to stay safe we choose min(1. / lambda_max, 1E-2)
            double omega;
            double eigen_value;
            Eigen::SparseMatrix <double> Iteration_Matrix = Diagonal * this->solver.Matrix_main.transpose() * this->solver.Matrix_main;
            Eigen::VectorXd eigen_vec_prev = Eigen::VectorXd::Ones(this->statistic.num_variable);
            eigen_vec_prev /= (eigen_vec_prev.array() * eigen_vec_prev.array()).sum();
            Eigen::VectorXd eigen_vec = Eigen::VectorXd::Ones(this->statistic.num_variable);
            while(true){
                eigen_vec = Iteration_Matrix * eigen_vec;
                eigen_value = (eigen_vec.array() * (Iteration_Matrix * eigen_vec).array()).sum();
                eigen_value /= (eigen_vec.array() * eigen_vec.array()).sum();
                eigen_vec /= eigen_value;

                if((eigen_vec - eigen_vec_prev).norm() >= 1E-12){
                    eigen_vec_prev = eigen_vec;
                }
                else{
                    break;
                }
            }
//            omega = std::min(1. / abs(eigen_value), 1E-3);
            omega = 1.;

            // Initialization
            double rho = 1E-4;
            this->solver.sol.prime.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.variables.now = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.prime.price_margin = Eigen::VectorXd::Zero(this->statistic.num_variable);
            this->solver.sol.dual.variables.prev = Eigen::VectorXd::Zero(this->statistic.num_constraint);
            this->solver.sol.dual.variables.now = Eigen::VectorXd::Zero(this->statistic.num_constraint);

            // Main loop
            // f(x) + .5 * rho * ||A %*% x + u||^2
            // f(x) + .5 * rho * ||a_i %*% x_i + (A_!i %*% x_!i + u)||^2
            // f(x) + .5 * rho * (a_i^2 x_i^2 + 2 * (b * a_i) %*% x_i + ...)
            // KKT: p_i(x) + rho * (a_i^2 x_i + (b * a_i)) = 0
            // or p_i(x) = -rho * (a_i^2 x_i + (b * a_i))
            // or x_i = (-p_i(x) / rho - b * a_i) / a_i^2
            int loop = 0;
            while(true){
                // Update prime variables and associate prices
                while(true){
//                    Eigen::VectorXd constant = this->solver.Matrix_main * this->solver.sol.prime.variables.now + this->solver.sol.dual.variables.now;
                    for(int var_iter = 0; var_iter < this->statistic.num_variable; ++ var_iter){
                        Eigen::VectorXd constant = this->solver.Matrix_main * this->solver.sol.prime.variables.now + this->solver.sol.dual.variables.now;
                        Eigen::VectorXd constant_temp = constant - this->solver.sol.prime.variables.now(var_iter) * this->solver.Matrix_main.col(var_iter);
                        double inner_prod = (constant_temp.transpose() * this->solver.Matrix_main.col(var_iter)).sum();

                        Eigen::Vector2i gap_price_ID(0, this->obj.cost_funcs[var_iter].moc.price.size() - 1);
                        int mid_price_ID = gap_price_ID.sum() / 2.;
                        while(gap_price_ID(1) - gap_price_ID(0) > 1){
                            Eigen::Vector2d gap_rent;
                            gap_rent(0) = this->obj.cost_funcs[var_iter].moc.price(gap_price_ID(0));
                            gap_rent(0) += rho * (sensitivity_x(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(0)) + inner_prod);
                            gap_rent(1) = this->obj.cost_funcs[var_iter].moc.price(gap_price_ID(1));
                            gap_rent(1) += rho * (sensitivity_x(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(gap_price_ID(1)) + inner_prod);

                            double mid_rent;
                            mid_rent = this->obj.cost_funcs[var_iter].moc.price(mid_price_ID);
                            mid_rent += rho * (sensitivity_x(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID) + inner_prod);

                            // Bisection root-search method
                            if(mid_rent * gap_rent(0) <= 0){
                                gap_price_ID(1) = mid_price_ID;
                            }
                            else{
                                gap_price_ID(0) = mid_price_ID;
                            }

                            mid_price_ID = gap_price_ID.sum() / 2.;
                        }

                        // Check binding term
                        if(gap_price_ID(0) % 2 == 1){
                            // price fixed case
                            this->solver.sol.prime.price_margin(var_iter) = this->obj.cost_funcs[var_iter].moc.price(mid_price_ID);
                            this->solver.sol.prime.variables.now(var_iter) = -this->solver.sol.prime.price_margin(var_iter) / rho;
                            this->solver.sol.prime.variables.now(var_iter) -= inner_prod;
                            this->solver.sol.prime.variables.now(var_iter) /= sensitivity_x(var_iter);
                        }
                        else{
                            // quantity fixed case
                            this->solver.sol.prime.variables.now(var_iter) = this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID);
                            this->solver.sol.prime.price_margin(var_iter) = -rho * (sensitivity_x(var_iter) * this->obj.cost_funcs[var_iter].moc.quantity(mid_price_ID) + inner_prod);
                        }
                    }

                    // Stabilization weight for x update
                    this->solver.sol.prime.variables.now = omega * this->solver.sol.prime.variables.now + (1. - omega) * this->solver.sol.prime.variables.prev;
                    this->solver.sol.prime.variables.prev = this->solver.sol.prime.variables.now;

                    // Check convergence for subproblem
                    double tol_sub = 10. / (loop + 1);
                    tol_sub = std::max(tol_sub, 1E-2 * tol_dual);
                    this->solver.sol.dual.error = this->solver.Matrix_main * this->solver.sol.prime.variables.now + this->solver.sol.dual.variables.now;
                    this->solver.sol.dual.error = rho * this->solver.Matrix_main.transpose() * this->solver.sol.dual.error;
                    this->solver.sol.dual.error += this->solver.sol.prime.price_margin;
                    if((this->solver.sol.dual.error.array() * this->solver.sol.dual.error.array()).maxCoeff() < tol_sub){
                        break;
                    }
                }

                // Update dual variables
                this->solver.sol.dual.variables.now += this->solver.Matrix_main * this->solver.sol.prime.variables.now;

                // Check convergence for whole problem
                this->solver.sol.prime.error = this->solver.Matrix_main * this->solver.sol.prime.variables.now;
                this->solver.sol.dual.error = this->solver.sol.prime.price_margin;
                this->solver.sol.dual.error += rho * this->solver.Matrix_main.transpose() * this->solver.sol.dual.variables.now;
                double prime_error_norm = this->solver.sol.prime.error.norm();
                prime_error_norm /= this->solver.sol.prime.error.size();
                double dual_error_norm = this->solver.sol.dual.error.norm();
                dual_error_norm /= this->solver.sol.dual.error.size();
                if(prime_error_norm < tol_prime && dual_error_norm < tol_dual){
                    break;
                }

                // Print progress
                int sys_show = (int) 6. - log(this->statistic.num_variable) / log(10.);
                sys_show = std::max(sys_show, 0);
                if(loop % (int) pow(10., sys_show) == 0){
                    std::cout << "Loop:\t" << loop << "\n";
                    std::cout << "Prime Error:\t" << prime_error_norm << "\n";
                    std::cout << "Dual Error:\t" << dual_error_norm << "\n";
                    std::cout << "\n\n";
                }

                // Update rho according to current prime and dual errors
                if(loop < 1E6){
                    rho += 1E-4;
                    rho = std::min(100., rho);
                }
                loop += 1;
            }
            std::cout << "Total loop:\t" << loop << "\n";
            std::cout << "Prime Error:\t" << this->solver.sol.prime.error.norm() / this->solver.sol.prime.error.size() << "\n";
            std::cout << "Dual Error:\t" << this->solver.sol.dual.error.norm() / this->solver.sol.dual.error.size() << "\n";
            std::cout << "Solution:\n";
            std::cout << this->solver.sol.prime.variables.now.segment(this->statistic.num_node, this->statistic.num_node).transpose() << "\n";
        }
    };

    // Functions
    void radial_line_problem_set(opf_struct&, int, int, std::complex<double>, double, double, double);
}




