// Functions for ADMM Algorithm
#include "header.h"

namespace ADMM{
    void radial_line_problem_split_set(opf_structs &opfs, int num_node, int num_price, std::complex<double> y_l, double theta_limit, double current_limit, double total_load){
        // Set full system statistic
        opfs.statistic.num_node = num_node;
        opfs.statistic.num_line = num_node - 1;
        opfs.opf_sub = std::vector <opf_struct> (3);
        int num_line_total = opfs.statistic.num_line;

        // Set first level of network
        {
            opfs.opf_sub[0].network.network_map.current_level = 0;
            opfs.opf_sub[0].statistic.num_node = 2;
            opfs.opf_sub[0].statistic.num_line = opfs.opf_sub[0].statistic.num_node - 1;
            int num_line = opfs.opf_sub[0].statistic.num_line;
            // Set network information
            opfs.opf_sub[0].network.line_conductance = y_l * Eigen::VectorXcd::Ones(num_line);
            opfs.opf_sub[0].network.line_conductance *= num_line_total;
            opfs.opf_sub[0].network.topology.reserve(num_line);
            opfs.opf_sub[0].network.topology.push_back(Eigen::Vector2i(0, 1));
        }

        // Set second level of network
        for(int network_iter = 1; network_iter < 3; ++ network_iter){
            bool left_vol_fix = network_iter == 2;
            bool right_vol_fix = network_iter == 1;

            // Set systems statistic
            opfs.opf_sub[network_iter].network.network_map.current_level = 1;
            opfs.opf_sub[network_iter].statistic.num_node = opfs.statistic.num_node / 2;
            if(network_iter == 1){
                opfs.opf_sub[network_iter].statistic.num_node += (opfs.statistic.num_node != (opfs.statistic.num_node / 2) * 2);
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

            // Connection nodes
            for(int node_iter = 0; node_iter < num_node; ++ node_iter){
                if(network_iter == 1 && node_iter == num_node){
                    opfs.opf_sub[network_iter].network.network_map.connection_nodes[node_iter] = 0;
                }
                else if(network_iter == 2 && node_iter == 0){
                    opfs.opf_sub[network_iter].network.network_map.connection_nodes[node_iter] = 1;
                }
                else{
                    opfs.opf_sub[network_iter].network.network_map.connection_nodes[node_iter] = -1;
                }
            }

            // Green functions
            int num_free_node = num_node - 1;
            opfs.opf_sub[network_iter].solver.constraint.green_func.voltage = Eigen::MatrixXd::Zero(num_free_node, num_free_node);
            opfs.opf_sub[network_iter].solver.constraint.green_func.current = Eigen::MatrixXd::Zero(num_line, num_free_node);
            Eigen::VectorXd z_cum(num_node);
            z_cum(0) = 0.;
            for(int node_iter = 1; node_iter < num_node; ++ node_iter){
                z_cum(node_iter) = z_cum(node_iter - 1) + 1. / y_l.imag();
            }

            if(network_iter < 2){
                for(int node_iter_1 = 0; node_iter_1 < num_free_node; ++ node_iter_1){
                    double z_L = z_cum(node_iter_1);
                    double z_R = z_cum.maxCoeff() - z_L;
                    double dV = z_R;

                    for(int node_iter_2 = 0; node_iter_2 < num_free_node; ++ node_iter_2){
                        if(node_iter_2 <= node_iter_1){
                            opfs.opf_sub[network_iter].solver.constraint.green_func.voltage(node_iter_2, node_iter_1) = dV;
                        }
                        else{
                            double z_L_temp = z_cum(node_iter_2);
                            double z_R_temp = z_cum.maxCoeff() - z_L_temp;

                            opfs.opf_sub[network_iter].solver.constraint.green_func.voltage(node_iter_2, node_iter_1) = dV * z_R_temp / z_L;
                        }
                    }

                    for(int line_iter = node_iter_1 + 1; line_iter < num_line; ++ line_iter){
                        opfs.opf_sub[network_iter].solver.constraint.green_func.current(line_iter, node_iter_1) = 1;
                    }
                }
            }
            else{
                for(int node_iter_1 = 0; node_iter_1 < num_free_node; ++ node_iter_1){
                    double z_L = z_cum(node_iter_1);
                    double dV = z_L;

                    for(int node_iter_2 = 0; node_iter_2 < num_free_node; ++ node_iter_2){
                        if(node_iter_2 <= node_iter_1){
                            double z_L_temp = z_cum(node_iter_2 + 1);

                            opfs.opf_sub[network_iter].solver.constraint.green_func.voltage(node_iter_2, node_iter_1) = dV * z_L_temp / z_L;
                        }
                        else{
                            opfs.opf_sub[network_iter].solver.constraint.green_func.voltage(node_iter_2, node_iter_1) = dV;
                        }
                    }

                    for(int line_iter = 0; line_iter < node_iter_1; ++ line_iter){
                        opfs.opf_sub[network_iter].solver.constraint.green_func.current(line_iter, node_iter_1) = -1;
                    }
                }
            }

            // Set cost coeff
            opfs.opf_sub[network_iter].solver.constraint.cost_coeff.voltage = 0.;
            opfs.opf_sub[network_iter].solver.constraint.cost_coeff.current = 0.;
            for(int node_iter_1 = 0; node_iter_1 < num_free_node; ++ node_iter_1){
                // Exchange with boundary
                for(int node_iter_2 = 0; node_iter_2 < num_free_node; ++ node_iter_2){
                    double value_temp = 1. / opfs.opf_sub[network_iter].solver.constraint.green_func.voltage(node_iter_2, node_iter_1);
                    value_temp = abs(value_temp);
                    opfs.opf_sub[network_iter].solver.constraint.cost_coeff.voltage = std::max(opfs.opf_sub[network_iter].solver.constraint.cost_coeff.voltage, value_temp);
                }

                for(int line_iter = 0; line_iter < num_line; ++ line_iter){
                    double value_temp;
                    if(opfs.opf_sub[network_iter].solver.constraint.green_func.current(line_iter, node_iter_1) != 0.){
                        value_temp = 1. / opfs.opf_sub[network_iter].solver.constraint.green_func.current(line_iter, node_iter_1);
                        value_temp = abs(value_temp);
                    }
                    else{
                        value_temp = 0.;
                    }
                    opfs.opf_sub[network_iter].solver.constraint.cost_coeff.current = std::max(opfs.opf_sub[network_iter].solver.constraint.cost_coeff.current, value_temp);
                }

                // Exchange between nodes
                for(int node_iter_2 = 0; node_iter_2 < num_free_node; ++ node_iter_2){
                    Eigen::VectorXd green_V_temp = opfs.opf_sub[network_iter].solver.constraint.green_func.voltage.col(node_iter_1);
                    green_V_temp -= opfs.opf_sub[network_iter].solver.constraint.green_func.voltage.col(node_iter_2);
                    Eigen::VectorXd green_I_temp = opfs.opf_sub[network_iter].solver.constraint.green_func.current.col(node_iter_1);
                    green_I_temp -= opfs.opf_sub[network_iter].solver.constraint.green_func.current.col(node_iter_2);

                    for(int node_iter_3 = 0; node_iter_3 < num_free_node; ++ node_iter_3){
                        double value_temp;
                        if(green_V_temp(node_iter_3) != 0.){
                            value_temp = 1. / green_V_temp(node_iter_3);
                            value_temp = abs(value_temp);
                        }
                        else{
                            value_temp = 0.;
                        }
                        opfs.opf_sub[network_iter].solver.constraint.cost_coeff.voltage = std::max(opfs.opf_sub[network_iter].solver.constraint.cost_coeff.voltage, value_temp);
                    }

                    for(int line_iter = 0; line_iter < num_line; ++ line_iter){
                        double value_temp;
                        if(green_I_temp(line_iter) != 0.){
                            value_temp = 1. / green_I_temp(line_iter);
                            value_temp = abs(value_temp);
                        }
                        else{
                            value_temp = 0.;
                        }
                        opfs.opf_sub[network_iter].solver.constraint.cost_coeff.current = std::max(opfs.opf_sub[network_iter].solver.constraint.cost_coeff.current, value_temp);
                    }
                }
            }

            // Set solver
            opfs.opf_sub[network_iter].transformation_set(0);
            opfs.opf_sub[network_iter].DC_Matrix_main_set(left_vol_fix, right_vol_fix);

            // Set cost functions
            opfs.opf_sub[network_iter].obj.cost_funcs = std::vector <opf_struct::obj_struct::cost_func_struct> (opfs.opf_sub[network_iter] .statistic.num_variable);
            opfs.opf_sub[network_iter].obj.price_range << -500., 3000.;
            double price_gap = opfs.opf_sub[network_iter].obj.price_range(1) - opfs.opf_sub[network_iter].obj.price_range(0);
            double quantity_inflex = 1000.;

            // Phase angle boundaries
            int num_voltage = num_node - left_vol_fix - right_vol_fix;
            for(int node_iter = 0; node_iter < num_voltage; ++ node_iter){
                int var_ID = node_iter;

                // Set bid functions for suuply
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.price = Eigen::VectorXd(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.quantity = Eigen::VectorXd::Zero(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.price << -std::numeric_limits<double>::infinity(), 0., price_gap * opfs.opf_sub[network_iter].solver.constraint.cost_coeff.voltage, std::numeric_limits<double>::infinity();
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.quantity << 0., theta_limit, quantity_inflex, 0.;

                // Set bid functions for demand
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.price = Eigen::VectorXd(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.quantity = Eigen::VectorXd(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.price << -std::numeric_limits<double>::infinity(), -price_gap * opfs.opf_sub[network_iter].solver.constraint.cost_coeff.voltage, 0., std::numeric_limits<double>::infinity();
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.quantity << 0., quantity_inflex, theta_limit, 0.;

                // Set merit order curve for residual load
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].moc_set();
            }

            // Line current boundaries
            for(int line_iter = 0; line_iter < num_line; ++ line_iter){
                int var_ID = num_voltage + num_node + line_iter;

                // Set bid functions for suuply
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.price = Eigen::VectorXd(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.quantity = Eigen::VectorXd::Zero(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.price << -std::numeric_limits<double>::infinity(), 0., price_gap * opfs.opf_sub[network_iter].solver.constraint.cost_coeff.current, std::numeric_limits<double>::infinity();
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].supply.quantity << 0., theta_limit, quantity_inflex, 0.;

                // Set bid functions for demand
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.price = Eigen::VectorXd(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.quantity = Eigen::VectorXd(4);
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.price << -std::numeric_limits<double>::infinity(), -price_gap * opfs.opf_sub[network_iter].solver.constraint.cost_coeff.current, 0., std::numeric_limits<double>::infinity();
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].demand.quantity << 0., quantity_inflex, theta_limit, 0.;

                // Set merit order curve for residual load
                opfs.opf_sub[network_iter].obj.cost_funcs[var_ID].moc_set();
            }

            // Power source / sink cost functions
            // First node is source node
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.price = Eigen::VectorXd(3);
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.quantity = Eigen::VectorXd(3);
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.price << -std::numeric_limits<double>::infinity(), 0., std::numeric_limits<double>::infinity();
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].supply.quantity << 0., total_load, 0.;
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].demand.price = Eigen::VectorXd(2);
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].demand.quantity = Eigen::VectorXd::Zero(2);
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].demand.price << -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();
            opfs.opf_sub[network_iter].obj.cost_funcs[num_node].moc_set();

            // Other nodes are sink node
            double node_load = total_load / (num_node - 1) / num_price;
            for(int node_iter = 1; node_iter < num_node; ++ node_iter){
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
        }
    }

    // Problem initialization for test simple radial power line
    void radial_line_problem_set(opf_struct &opf, int num_node, int num_price, std::complex<double> y_l, double theta_limit, double current_limit, double total_load){
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

        // Set solver
        opf.transformation_set(0);
        opf.DC_Matrix_main_set();

        // Set cost functions
        opf.obj.cost_funcs = std::vector <opf_struct::obj_struct::cost_func_struct> (opf.statistic.num_variable);

        // Phase angle boundaries
        for(int node_iter = 0; node_iter < num_node; ++ node_iter){
            int var_ID = node_iter;
            opf_struct::obj_struct::cost_func_struct cost_func;
            cost_func.moc.price = Eigen::VectorXd(4);
            cost_func.moc.quantity = Eigen::VectorXd(4);
            cost_func.moc.obj = Eigen::VectorXd::Zero(4);

            cost_func.moc.price << -std::numeric_limits<double>::infinity(), 0., 0., std::numeric_limits<double>::infinity();
            cost_func.moc.quantity << -theta_limit, -theta_limit, theta_limit, theta_limit;

            opf.obj.cost_funcs[var_ID] = cost_func;
        }

        // Line current boundaries
        for(int line_iter = 0; line_iter < num_line; ++ line_iter){
            int var_ID = 2 * num_node + line_iter;
            opf_struct::obj_struct::cost_func_struct cost_func;
            cost_func.moc.price = Eigen::VectorXd(4);
            cost_func.moc.quantity = Eigen::VectorXd(4);
            cost_func.moc.obj = Eigen::VectorXd::Zero(4);

            cost_func.moc.price << -std::numeric_limits<double>::infinity(), 0., 0., std::numeric_limits<double>::infinity();
            cost_func.moc.quantity << -current_limit, -current_limit, current_limit, current_limit;

            opf.obj.cost_funcs[var_ID] = cost_func;
        }

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

        // Other nodes are sink node
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
    }
}
