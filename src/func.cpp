// Functions for ADMM Algorithm
#include "header.h"

namespace ADMM{
    void radial_line_problem_split_set(opf_structs &opfs, int num_node, int num_price, std::complex<double> y_l, double theta_limit, double current_limit, double total_load, double penalty_price_voltage){
        // Set full system statistic
        opfs.statistic.num_node = num_node;
        opfs.statistic.num_line = num_node - 1;
        opfs.opf_sub = std::vector <opf_struct> (2);
        int num_line_total = opfs.statistic.num_line;

        // Initialize subnetworks
        // Use hierarchy topology information later if n-split actually is faster
        for(int network_iter = 0; network_iter < opfs.opf_sub.size(); ++ network_iter){
            // Set systems statistic
            opfs.opf_sub[network_iter].statistic.num_node = opfs.statistic.num_node / opfs.opf_sub.size();
            if(network_iter == opfs.opf_sub.size()){
                opfs.opf_sub[network_iter].statistic.num_node += opfs.statistic.num_node % opfs.opf_sub.size();
            }
            opfs.opf_sub[network_iter].statistic.num_line = opfs.opf_sub[network_iter].statistic.num_node - 1;
            opfs.opf_sub[network_iter].problem_size_parmeters_set();
            int num_node = opfs.opf_sub[network_iter].statistic.num_node;
            int num_line = opfs.opf_sub[network_iter].statistic.num_line;

            // Set network information
            opfs.opf_sub[network_iter].network.line_conductance = y_l * Eigen::VectorXcd::Ones(num_line);
            opfs.opf_sub[network_iter].network.line_conductance *= num_line_total;
            opfs.opf_sub[network_iter].network.topology.reserve(num_line);
            for(int line_iter = 0; line_iter < num_line; ++ line_iter){
                opfs.opf_sub[network_iter].network.topology.push_back(Eigen::Vector2i(line_iter, line_iter + 1));
            }

            // Set cost functions and merit order curves
            opfs.opf_sub[network_iter].moc_initialize(theta_limit, current_limit, penalty_price_voltage);

            // Power source / sink cost functions
            // First node is source node
            if(network_iter == 0){
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.price = Eigen::VectorXd(3);
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.quantity = Eigen::VectorXd(3);
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.price << -std::numeric_limits<double>::infinity(), 0., std::numeric_limits<double>::infinity();
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.quantity << 0., total_load, 0.;
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].demand.price = Eigen::VectorXd(2);
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].demand.quantity = Eigen::VectorXd::Zero(2);
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].demand.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
                opfs.opf_sub[network_iter].obj.cost_funcs[num_node].moc_set();
            }

            // Other nodes are sink nodes
            double node_load = total_load / (num_node - 1) / num_price;
            for(int node_iter = (network_iter == 0); node_iter < num_node; ++ node_iter){
                int var_ID = num_node + node_iter;

                // Set bid functions for suuply
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.price = Eigen::VectorXd(2);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.quantity = Eigen::VectorXd::Zero(2);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();

                // Set bid functions for demand
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.price = Eigen::VectorXd(num_price + 2);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.quantity = Eigen::VectorXd(num_price + 2);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.price << -std::numeric_limits<double>::infinity(), Eigen::VectorXd::LinSpaced(num_price, 1., (double) num_price), std::numeric_limits<double>::infinity();
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.quantity << 0., node_load * Eigen::VectorXd::Ones(num_price), 0.;

                // Set merit order curve for residual load
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].moc_set();
            }

            // Set solver
            opfs.opf_sub[network_iter].transformation_set(0);
            opfs.opf_sub[network_iter].DC_Matrix_main_set();
        }
    }

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

        // Set cost functions and merit order curves
        opf.moc_initialize(theta_limit, current_limit, penalty_price_voltage);

        // Power source / sink cost functions
        // First node is source node
        opf.obj.cost_funcs[num_node].supply.price = Eigen::VectorXd(3);
        opf.obj.cost_funcs[num_node].supply.quantity = Eigen::VectorXd(3);
        opf.obj.cost_funcs[num_node].supply.price << -std::numeric_limits<double>::infinity(), 0., std::numeric_limits<double>::infinity();
        opf.obj.cost_funcs[num_node].supply.quantity << 0., total_load, 0.;
        opf.obj.cost_funcs[num_node].demand.price = Eigen::VectorXd(2);
        opf.obj.cost_funcs[num_node].demand.quantity = Eigen::VectorXd::Zero(2);
        opf.obj.cost_funcs[num_node].demand.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
        opf.obj.cost_funcs[num_node].moc_set();

        // Other nodes are sink nodes
        double node_load = total_load / (num_node - 1) / num_price;
        for(int node_iter = 1; node_iter < num_node; ++ node_iter){
            int var_ID = num_node + node_iter;

            // Set bid functions for suuply
            opf.obj.cost_funcs[var_ID].supply.price = Eigen::VectorXd(2);
            opf.obj.cost_funcs[var_ID].supply.quantity = Eigen::VectorXd::Zero(2);
            opf.obj.cost_funcs[var_ID].supply.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();

            // Set bid functions for demand
            opf.obj.cost_funcs[var_ID].demand.price = Eigen::VectorXd(num_price + 2);
            opf.obj.cost_funcs[var_ID].demand.quantity = Eigen::VectorXd(num_price + 2);
            opf.obj.cost_funcs[var_ID].demand.price << -std::numeric_limits<double>::infinity(), Eigen::VectorXd::LinSpaced(num_price, 1., (double) num_price), std::numeric_limits<double>::infinity();
            opf.obj.cost_funcs[var_ID].demand.quantity << 0., node_load * Eigen::VectorXd::Ones(num_price), 0.;

            // Set merit order curve for residual load
            opf.obj.cost_funcs[var_ID].moc_set();
        }

        // Set solver
        opf.transformation_set(0);
        opf.DC_Matrix_main_set();
    }
}
