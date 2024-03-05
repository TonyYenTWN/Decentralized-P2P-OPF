// Functions for ADMM Algorithm
#include "header.h"

namespace ADMM{
    void radial_line_problem_set(opf_struct &opf, int num_node, int num_price, std::complex<double> y_l, double theta_limit, double current_limit, double total_load, double penalty_price_voltage){
        // Set systems statistic
        opf.statistic.num_node = num_node;
        opf.statistic.num_line = num_node - 1;
        int num_line = opf.statistic.num_line;
        opf.problem_size_parmeters_set();

        // Set network information
        opf.network.line_conductance = y_l * Eigen::VectorXcd::Ones(num_line);
        opf.network.line_conductance *= num_line;
        opf.network.topology.reserve(num_line);
        for(int line_iter = 0; line_iter < num_line; ++ line_iter){
            opf.network.topology.push_back(Eigen::Vector2i(line_iter, line_iter + 1));
        }

        // Power source / sink cost functions
        // First node is source node
        opf.moc_initialize(current_limit);
        opf.obj.cost_funcs[0].supply.price = Eigen::VectorXd(3);
        opf.obj.cost_funcs[0].supply.quantity = Eigen::VectorXd(3);
        opf.obj.cost_funcs[0].supply.price << -std::numeric_limits<double>::infinity(), 0., std::numeric_limits<double>::infinity();
        opf.obj.cost_funcs[0].supply.quantity << 0., total_load, 0.;
        opf.obj.cost_funcs[0].demand.price = Eigen::VectorXd(2);
        opf.obj.cost_funcs[0].demand.quantity = Eigen::VectorXd::Zero(2);
        opf.obj.cost_funcs[0].demand.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
        opf.obj.cost_funcs[0].moc_set();

        // Other nodes are sink nodes
        double node_load = total_load / (num_node - 1) / num_price;
        for(int node_iter = 1; node_iter < num_node; ++ node_iter){
            // Set bid functions for suuply
            opf.obj.cost_funcs[node_iter].supply.price = Eigen::VectorXd(2);
            opf.obj.cost_funcs[node_iter].supply.quantity = Eigen::VectorXd::Zero(2);
            opf.obj.cost_funcs[node_iter].supply.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();

            // Set bid functions for demand
            opf.obj.cost_funcs[node_iter].demand.price = Eigen::VectorXd(num_price + 2);
            opf.obj.cost_funcs[node_iter].demand.quantity = Eigen::VectorXd(num_price + 2);
            opf.obj.cost_funcs[node_iter].demand.price << -std::numeric_limits<double>::infinity(), Eigen::VectorXd::LinSpaced(num_price, 1., (double) num_price), std::numeric_limits<double>::infinity();
            opf.obj.cost_funcs[node_iter].demand.quantity << 0., node_load * Eigen::VectorXd::Ones(num_price), 0.;

            // Set merit order curve for residual load
            opf.obj.cost_funcs[node_iter].moc_set();
        }

        // Set solver
        opf.DC_Matrix_main_set();
    }
}
