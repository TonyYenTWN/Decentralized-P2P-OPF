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
        struct bid_struct{
            Eigen::MatrixXd demand;
            Eigen::MatrixXd supply;

            Eigen::MatrixXd moc_rl_set(){
//                Eigen::MatrixXd moc_rl
            }
        };
        struct obj_struct{
            double theta_mul;
            std::vector <bid_struct> bid_node;
            std::vector <Eigen::MatrixXd> cost_funcs;
        };
        obj_struct obj;

        // Solver
        struct solver_struct{
            // Main matrix
            Eigen::SparseMatrix <double> Matrix_main;
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




