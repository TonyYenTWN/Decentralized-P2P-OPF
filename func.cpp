// Functions for ADMM Algorithm
#include "header.h"

namespace ADMM{
    void radial_line_problem_set(opf_struct &opf, int num_node, int num_price, std::complex<double> y_l, double theta_limit, double current_limit, double total_load){
        // Set systems statistic
        opf.statistic.num_node = num_node;
        opf.statistic.num_line = num_node - 1;
        int num_line = opf.statistic.num_line;
        opf.problem_size_parmeters_set();

        // Set network information
        opf.network.line_conductance = y_l * Eigen::VectorXcd::Ones(num_line);
        opf.network.line_conductance /= num_line;
        opf.network.topology.reserve(num_line);
        for(int line_iter = 0; line_iter < num_line; ++ line_iter){
            opf.network.topology.push_back(Eigen::Vector2i(line_iter, line_iter + 1));
        }

        // Set cost functions
        opf.theta_mul_calc();
        opf.obj.cost_funcs = std::vector <Eigen::MatrixXd> (opf.statistic.num_variable);

        // Phase angle boundaries
        for(int node_iter = 0; node_iter < num_node; ++ node_iter){
            int var_ID = node_iter;
            Eigen::MatrixXd cost_func(4, 2);
            cost_func.row(0) << -std::numeric_limits<double>::infinity(), -theta_limit * opf.obj.theta_mul;
            cost_func.row(1) << 0, -theta_limit * opf.obj.theta_mul;
            cost_func.row(2) << 0, theta_limit * opf.obj.theta_mul;
            cost_func.row(3) << std::numeric_limits<double>::infinity(), theta_limit * opf.obj.theta_mul;
            opf.obj.cost_funcs[var_ID] = cost_func;
        }

        // Line current boundaries
        for(int line_iter = 0; line_iter < num_line; ++ line_iter){
            int var_ID = 2 * num_node + line_iter;
            Eigen::MatrixXd cost_func(4, 2);
            cost_func.row(0) << -std::numeric_limits<double>::infinity(), -current_limit;
            cost_func.row(1) << 0, -theta_limit * opf.obj.theta_mul;
            cost_func.row(2) << 0, theta_limit * opf.obj.theta_mul;
            cost_func.row(3) << std::numeric_limits<double>::infinity(), current_limit;
            opf.obj.cost_funcs[var_ID] = cost_func;
        }

        // Power source / sink cost functions
        opf.obj.bid_node = std::vector <opf_struct::obj_struct::bid_struct> (num_node);

        // First node is source node
        opf.obj.bid_node[0].supply.price = Eigen::VectorXd(3);
        opf.obj.bid_node[0].supply.quantity = Eigen::VectorXd(3);
        opf.obj.bid_node[0].supply.price << -std::numeric_limits<double>::infinity(), 0., std::numeric_limits<double>::infinity();
        opf.obj.bid_node[0].supply.quantity << 0., total_load, 0.;
        opf.obj.bid_node[0].demand.price = Eigen::VectorXd(2);
        opf.obj.bid_node[0].demand.quantity = Eigen::VectorXd::Zero(2);
        opf.obj.bid_node[0].demand.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
        opf.obj.bid_node[0].moc_rl_set();
        std::cout << "\n";

        // Other nodes are sink node
        double node_load = total_load / (num_node - 1) / num_price;
        for(int node_iter = 1; node_iter < num_node; ++ node_iter){
            int var_ID = num_node + node_iter;

            // Set bid functions for suuply
            opf.obj.bid_node[node_iter].supply.price = Eigen::VectorXd(2);
            opf.obj.bid_node[node_iter].supply.quantity = Eigen::VectorXd::Zero(2);
            opf.obj.bid_node[node_iter].supply.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();

            // Set bid functions for demand
            opf.obj.bid_node[node_iter].demand.price = Eigen::VectorXd(num_price + 2);
            opf.obj.bid_node[node_iter].demand.quantity = Eigen::VectorXd(num_price + 2);
            opf.obj.bid_node[node_iter].demand.price << -std::numeric_limits<double>::infinity(), Eigen::VectorXd::LinSpaced(num_price, 1., (double) num_price), std::numeric_limits<double>::infinity();
            opf.obj.bid_node[node_iter].demand.quantity << 0., node_load * Eigen::VectorXd::Ones(num_price), 0.;

            // Set merit order curve for residual load
            opf.obj.bid_node[node_iter].moc_rl_set();
            std::cout << "\n";
        }

        // Set solver
        opf.Matrix_main_set();
    }


}
