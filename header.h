// header for ADMM

#pragma once
// STL
#include <complex>
#include <vector>
#include <iostream>

// Additional libraries
#include <boost/math/constants/constants.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Structures
namespace ADMM{
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
            struct bid_struct{
                // Substructures for bids
                struct bid_mini_struct{
                    Eigen::VectorXd price;
                    Eigen::VectorXd quantity;
                };
                bid_mini_struct demand;
                bid_mini_struct supply;

                // Substructure for merit order curve
                struct moc_struct{
                    Eigen::VectorXd price;
                    Eigen::VectorXd rl;
                    Eigen::VectorXd obj;
                };
                moc_struct moc_rl;

                void moc_rl_set(){
                    int num_row = 2 * (this->demand.price.size() + this->supply.price.size() - 3);
                    this->moc_rl.price = Eigen::VectorXd(num_row);
                    this->moc_rl.rl = Eigen::VectorXd(num_row);
                    this->moc_rl.obj = Eigen::VectorXd(num_row);

                    // Initialization for the loop
                    int price_demand_ID = 0;
                    int price_supply_ID = 0;
                    int price_rl_ID = 0;
                    double obj_value = -(this->demand.price.segment(1, this->demand.price.size() - 2).array() * this->demand.quantity.segment(1, this->demand.quantity.size() - 2).array()).sum();

                    // First item when price = -Inf
                    this->moc_rl.price(0) = -std::numeric_limits<double>::infinity();
                    this->moc_rl.rl(0) = -this->demand.quantity.sum();
                    this->moc_rl.obj(0) = obj_value;
                    std::cout << moc_rl.price(price_rl_ID) << "\t";
                    std::cout << moc_rl.rl(price_rl_ID) << "\t";
                    std::cout << moc_rl.obj(price_rl_ID) << "\t";
                    std::cout << "\n";

                    while(price_rl_ID < num_row - 2){
                        price_rl_ID += 1;

                        // Next demand marginal price lower than supply; add that demand quantity to residual load
                        if(this->demand.price(price_demand_ID + 1) <= this->supply.price(price_supply_ID + 1)){
                            price_demand_ID += price_demand_ID < this->demand.price.size() - 2;
                            this->moc_rl.price(price_rl_ID) = this->demand.price(price_demand_ID);
                            this->moc_rl.rl(price_rl_ID) = this->moc_rl.rl(price_rl_ID - 1);
                            this->moc_rl.obj(price_rl_ID, 2) = obj_value;
                            std::cout << moc_rl.price(price_rl_ID) << "\t";
                            std::cout << moc_rl.rl(price_rl_ID) << "\t";
                            std::cout << moc_rl.obj(price_rl_ID) << "\t";
                            std::cout << "\n";

                            price_rl_ID += 1;
                            obj_value += this->demand.price(price_demand_ID) * this->demand.quantity(price_demand_ID);
                            this->moc_rl.price(price_rl_ID) = this->demand.price(price_demand_ID);
                            this->moc_rl.rl(price_rl_ID) = this->moc_rl.rl(price_rl_ID - 1) + this->demand.quantity(price_demand_ID);
                            this->moc_rl.obj(price_rl_ID) = obj_value;
                            std::cout << moc_rl.price(price_rl_ID) << "\t";
                            std::cout << moc_rl.rl(price_rl_ID) << "\t";
                            std::cout << moc_rl.obj(price_rl_ID) << "\t";
                            std::cout << "\n";

                        }
                        // Next supply marginal price lower than demand; add that supply quantity to residual load
                        else{
                            price_supply_ID += price_supply_ID < this->supply.price.size() - 2;
                            this->moc_rl.price(price_rl_ID) = this->supply.price(price_supply_ID);
                            this->moc_rl.rl(price_rl_ID) = this->moc_rl.rl(price_rl_ID - 1);
                            this->moc_rl.obj(price_rl_ID, 2) = obj_value;
                            std::cout << moc_rl.price(price_rl_ID) << "\t";
                            std::cout << moc_rl.rl(price_rl_ID) << "\t";
                            std::cout << moc_rl.obj(price_rl_ID) << "\t";
                            std::cout << "\n";

                            price_rl_ID += 1;
                            obj_value += this->supply.price(price_supply_ID) * this->supply.quantity(price_supply_ID);
                            this->moc_rl.price(price_rl_ID) = this->supply.price(price_supply_ID);
                            this->moc_rl.rl(price_rl_ID) = this->moc_rl.rl(price_rl_ID - 1) + this->supply.quantity(price_supply_ID);
                            this->moc_rl.obj(price_rl_ID) = obj_value;
                            std::cout << moc_rl.price(price_rl_ID) << "\t";
                            std::cout << moc_rl.rl(price_rl_ID) << "\t";
                            std::cout << moc_rl.obj(price_rl_ID) << "\t";
                            std::cout << "\n";
                        }
                    }

                    // Last item when price = Inf
                    price_rl_ID += 1;
                    this->moc_rl.price(price_rl_ID) = std::numeric_limits<double>::infinity();
                    this->moc_rl.rl(price_rl_ID) = this->moc_rl.rl(price_rl_ID - 1);
                    this->moc_rl.obj(price_rl_ID) = obj_value;
                    std::cout << moc_rl.price(price_rl_ID) << "\t";
                    std::cout << moc_rl.rl(price_rl_ID) << "\t";
                    std::cout << moc_rl.obj(price_rl_ID) << "\t";
                    std::cout << "\n";
                }
            };

            double theta_mul;
            std::vector <bid_struct> bid_node;
            std::vector <Eigen::MatrixXd> cost_funcs;
        };
        obj_struct obj;

        // Solver
        struct solver_struct{
            // Main matrix
            Eigen::SparseMatrix <double> Matrix_main;

            // Solution
            struct sol_struct{
                struct value_struct{
                    Eigen::VectorXd variables;
                    double error;
                };

                value_struct prime;
                value_struct dual;
                double obj_value;
            };
        };
        solver_struct solver;

        // Functions
        void problem_size_parmeters_set(){
            this->statistic.num_variable = 2 * this->statistic.num_node + this->statistic.num_line;
            this->statistic.num_constraint = this->statistic.num_node + this->statistic.num_line;
        }

        void theta_mul_calc(){
            double log_y_l = 0;
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                log_y_l += log(std::abs(this->network.line_conductance(line_iter)));
            }
            this->obj.theta_mul = exp(log_y_l / this->statistic.num_line);
        }

        void Matrix_main_set(){
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
    };

    // Functions
    void radial_line_problem_set(opf_struct&, int, int, std::complex<double>, double, double, double);
}




