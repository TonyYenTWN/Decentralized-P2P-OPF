// header for ADMM

#pragma once

// STL
#include <complex>
#include <iostream>
#include <vector>

// Additional libraries
#include <Eigen/Dense>

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

        void rescale(Eigen::VectorXd scale){
             for(int row_iter = 0; row_iter < this->row.size(); ++ row_iter){
                for(int entry_iter = 0; entry_iter < this->row[row_iter].size(); ++ entry_iter){
                    int col_ID = this->row[row_iter][entry_iter].first;
                    this->row[row_iter][entry_iter].second *= scale(col_ID);
                }
             }

             for(int col_iter = 0; col_iter < this->col.size(); ++ col_iter){
                for(int entry_iter = 0; entry_iter < this->col[col_iter].size(); ++ entry_iter){
                    this->col[col_iter][entry_iter].second *= scale(col_iter);
                }
             }
        }

        SparseMatrix col_cross_prod(){
            SparseMatrix Mat;
            Mat.dim_set(this->col.size(), this->col.size());

            for(int row_iter = 0; row_iter < Mat.row.size(); ++ row_iter){
                for(int col_iter = 0; col_iter < Mat.col.size(); ++ col_iter){
                    int counter_0 = 0;
                    int counter_1 = 0;
                    double coeff = 0.;

                    bool loop_flag = counter_0 < this->col[row_iter].size();
                    loop_flag *= counter_1 < this->col[col_iter].size();
                    while(loop_flag){
                        // Add the multiplcation if entry is the same
                        if(this->col[row_iter][counter_0].first == this->col[col_iter][counter_1].first){
                            coeff += this->col[row_iter][counter_0].second * this->col[col_iter][counter_1].second;
                            counter_0 += 1;
                            counter_1 += 1;
                        }
                        else{
                            // Increase the smaller entry
                            if(this->col[row_iter][counter_0].first < this->col[col_iter][counter_1].first){
                                counter_0 += 1;
                            }
                            else{
                                counter_1 += 1;
                            }
                        }

                        // Update flag criteria
                        loop_flag = counter_0 < this->col[row_iter].size();
                        loop_flag *= counter_1 < this->col[col_iter].size();
                    }

                    if(coeff != 0.){
                        Mat.row[row_iter].push_back(std::pair <int, double> (col_iter, coeff));
                        Mat.col[col_iter].push_back(std::pair <int, double> (row_iter, coeff));
                    }
                }
            }

            return Mat;
        }

        Eigen::VectorXd col_vec_prod(Eigen::VectorXd vec){
            Eigen::VectorXd prod = Eigen::VectorXd::Zero(vec.size());

            for(int row_iter = 0; row_iter < this->row.size(); ++ row_iter){
                for(int entry_iter = 0; entry_iter < this->row[row_iter].size(); ++ entry_iter){
                    int col_ID = this->row[row_iter][entry_iter].first;
                    double coeff = this->row[row_iter][entry_iter].second;
                    prod(row_iter) += coeff * vec(col_ID);
                }
            }

            return prod;
        }
    };

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
            std::vector <Eigen::Vector2i> topology_line;
            std::vector <std::vector <Eigen::Vector2i>> topology_node;

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

                // Function to set merit order curve
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
                SparseMatrix Mat_main;
                SparseMatrix PSD_main;
                Eigen::VectorXd boundary;
            };
            constraint_struct constraint;

            // Infomation at each local (nodal) optimizer
            struct local_struct{
                struct ID_struct{
                    std::vector <int> variable;
                    std::vector <int> constraint;
                };
                ID_struct ID;

                struct information_struct{
                    struct relation_struct{
                        double self;
                        std::vector <std::pair <int, double>> neighbor;
                        std::vector <std::pair <int, double>> outside;
                    };
                    relation_struct prime;
                    relation_struct dual;
                };

                struct informations_struct{
                    std::vector <information_struct> variables;
                    std::vector <information_struct> constraints;
                };
                informations_struct infos;
            };
            std::vector <local_struct> locals;

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

        // Initialize merit order curves and set values for voltage and current
        void moc_initialize(double theta_limit, double current_limit){
            this->obj.cost_funcs = std::vector <obj_struct::cost_func_struct> (this->statistic.num_variable);

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

            // Line current boundaries
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int var_ID = 2 * this->statistic.num_node + line_iter;
                opf_struct::obj_struct::cost_func_struct cost_func;
                cost_func.moc.price = Eigen::VectorXd(4);
                cost_func.moc.quantity = Eigen::VectorXd(4);
                cost_func.moc.obj = Eigen::VectorXd::Zero(4);

                cost_func.moc.price << -std::numeric_limits<double>::infinity(), 0., 0., std::numeric_limits<double>::infinity();
                cost_func.moc.quantity << -current_limit, -current_limit, current_limit, current_limit;

                this->obj.cost_funcs[var_ID] = cost_func;
            }
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

        // Main and PSD matrix initialization for DC OPF
        void DC_Matrix_main_set(){
            // Initialize main matrix dimensions
            this->solver.constraint.Mat_main.dim_set(this->statistic.num_constraint, this->statistic.num_variable);

            // Node Balance Equation
            // S - t(NL) %*% I = 0
            for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                int row_ID = node_iter;

                {
                    int col_ID = this->statistic.num_node + node_iter;
                    double coeff = 1.;

                    this->solver.constraint.Mat_main.row[row_ID].push_back(std::pair <int, double> (col_ID, coeff));
                    this->solver.constraint.Mat_main.col[col_ID].push_back(std::pair <int, double> (row_ID, coeff));
                }

                for(int link_iter = 0; link_iter < this->network.topology_node[node_iter].size(); ++ link_iter){
                    int col_ID = 2 * this->statistic.num_node + this->network.topology_node[node_iter][link_iter](0);
                    double coeff = (double) -this->network.topology_node[node_iter][link_iter](1);

                    this->solver.constraint.Mat_main.row[row_ID].push_back(std::pair <int, double> (col_ID, coeff));
                    this->solver.constraint.Mat_main.col[col_ID].push_back(std::pair <int, double> (row_ID, coeff));
                }
            }

            // Line Current Equation
            // Y_l %*% NL %*% V - I = 0
            for(int line_iter = 0; line_iter < this->statistic.num_line; ++ line_iter){
                int row_ID = this->statistic.num_node + line_iter;
                double y_l = this->network.line_conductance(line_iter).imag();

                for(int end_iter = 0; end_iter < 2; ++ end_iter){
                    int col_ID = this->network.topology_line[line_iter](end_iter);
                    double coeff = y_l;
                    coeff *= 1 - 2 * end_iter;

                    this->solver.constraint.Mat_main.row[row_ID].push_back(std::pair <int, double> (col_ID, coeff));
                    this->solver.constraint.Mat_main.col[col_ID].push_back(std::pair <int, double> (row_ID, coeff));
                }

                {
                    int col_ID = 2 * this->statistic.num_node + line_iter;
                    double coeff = -1.;

                    this->solver.constraint.Mat_main.row[row_ID].push_back(std::pair <int, double> (col_ID, coeff));
                    this->solver.constraint.Mat_main.col[col_ID].push_back(std::pair <int, double> (row_ID, coeff));
                }
            }

            // Set boundary for equality constraints; transformation considered
            // Ax = c
            // A(mx' + x_0) = c
            // [A][m]x' = c - A * x_0
            this->solver.constraint.boundary = this->solver.constraint.Mat_main.col_vec_prod(-this->obj.transformation.shift);

            // Rescale Main Matrix
            this->solver.constraint.Mat_main.rescale(this->obj.transformation.scale);

            // Create PSD matrix
            this->solver.constraint.PSD_main = this->solver.constraint.Mat_main.col_cross_prod();
        }

        // Local information initialization
        void local_information_set(){
            this->solver.locals = std::vector <solver_struct::local_struct> (this->statistic.num_node);

            for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                int num_line_temp = this->network.topology_node[node_iter].size();
                int num_var_temp = 2 + num_line_temp;
                int num_constr_temp = 1 + num_line_temp;

                // Store the IDs
                int V_ID = node_iter;
                int S_ID = this->statistic.num_node + node_iter;
                this->solver.locals[node_iter].ID.variable.reserve(num_var_temp);
                this->solver.locals[node_iter].ID.constraint.reserve(num_constr_temp);
                this->solver.locals[node_iter].ID.variable.push_back(V_ID);
                this->solver.locals[node_iter].ID.constraint.push_back(V_ID);
                this->solver.locals[node_iter].ID.variable.push_back(S_ID);
                for(int link_iter = 0; link_iter < this->network.topology_node[node_iter].size(); ++ link_iter){
                    if(this->network.topology_node[node_iter][link_iter](1) != 1){
                        continue;
                    }

                    int I_ID = 2 * this->statistic.num_node + this->network.topology_node[node_iter][link_iter](0);
                    int constr_ID = this->statistic.num_node + this->network.topology_node[node_iter][link_iter](0);
                    this->solver.locals[node_iter].ID.variable.push_back(I_ID);
                    this->solver.locals[node_iter].ID.constraint.push_back(constr_ID);
                }

                // Store the information
                num_var_temp = this->solver.locals[node_iter].ID.variable.size();
                num_constr_temp = this->solver.locals[node_iter].ID.constraint.size();

                // Prime variables
                this->solver.locals[node_iter].infos.variables.reserve(num_var_temp);
                for(int var_iter = 0; var_iter < num_var_temp; ++ var_iter){
                    int var_ID = this->solver.locals[node_iter].ID.variable[var_iter];

                    // Check non-zero entries in PSD (for prime variables minimization)
                    this->solver.locals[node_iter].infos.variables[var_iter].prime.neighbor.reserve(num_var_temp);
                    this->solver.locals[node_iter].infos.variables[var_iter].prime.outside.reserve(this->statistic.num_variable);
                    for(int col_iter = 0; col_iter < this->solver.constraint.PSD_main.row[var_ID].size(); ++ col_iter){
                        auto pair_temp = this->solver.constraint.PSD_main.row[var_ID][col_iter];

                        bool break_flag = 0;
                        int var_temp_now = 0;
                        for(int var_temp_iter = var_temp_now; var_temp_iter < num_var_temp; ++ var_temp_iter){
                            int col_ID = pair_temp.first;
                            int var_ID_temp = this->solver.locals[node_iter].ID.variable[var_temp_iter];

                            if(col_ID == var_ID_temp){
                                if(col_ID == var_ID){
                                    double coeff = pair_temp.second;
                                    this->solver.locals[node_iter].infos.variables[var_iter].prime.self = coeff;
                                }
                                else{
                                    this->solver.locals[node_iter].infos.variables[var_iter].prime.neighbor.push_back(pair_temp);
                                }

                                var_temp_now = var_temp_iter + 1;

                                break_flag = 1;
                                break;
                            }
                        }

                        // No break occur in 2nd loop = outside
                        if(!break_flag){
                            this->solver.locals[node_iter].infos.variables[var_iter].prime.outside.push_back(pair_temp);
                        }
                    }

                    // Check non-zero entries in transpose of main (for prime variables minimization)
                    this->solver.locals[node_iter].infos.variables[var_iter].dual.neighbor.reserve(num_constr_temp);
                    this->solver.locals[node_iter].infos.variables[var_iter].dual.outside.reserve(this->statistic.num_constraint);
                    for(int row_iter = 0; row_iter < this->solver.constraint.Mat_main.col[var_ID].size(); ++ row_iter){
                        auto pair_temp = this->solver.constraint.Mat_main.col[var_ID][row_iter];

                        bool break_flag = 0;
                        int constr_temp_now = 0;
                        for(int constr_temp_iter = constr_temp_now; constr_temp_iter < num_constr_temp; ++ constr_temp_iter){
                            int row_ID = pair_temp.first;
                            int constr_ID_temp = this->solver.locals[node_iter].ID.constraint[constr_temp_iter];

                            if(row_ID == constr_ID_temp){
                                this->solver.locals[node_iter].infos.variables[var_iter].dual.neighbor.push_back(pair_temp);

                                constr_temp_now = constr_temp_iter + 1;

                                break_flag = 1;
                                break;
                            }
                        }

                        // No break occur in 2nd loop = outside
                        if(!break_flag){
                            this->solver.locals[node_iter].infos.variables[var_iter].dual.outside.push_back(pair_temp);
                        }
                    }
                }

                // Dual variables
                this->solver.locals[node_iter].infos.constraints.reserve(num_constr_temp);
                for(int constr_iter = 0; constr_iter < num_constr_temp; ++ constr_iter){
                    int constr_ID = this->solver.locals[node_iter].ID.constraint[constr_iter];

                    // Check non-zero entries in main (for dual variables maximization)
                    this->solver.locals[node_iter].infos.constraints[constr_iter].prime.neighbor.reserve(num_var_temp);
                    this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside.reserve(this->statistic.num_variable);
                    for(int col_iter = 0; col_iter < this->solver.constraint.Mat_main.row[constr_ID].size(); ++ col_iter){
                        auto pair_temp = this->solver.constraint.Mat_main.row[constr_ID][col_iter];

                        bool break_flag = 0;
                        int var_temp_now = 0;
                        for(int var_temp_iter = var_temp_now; var_temp_iter < num_var_temp; ++ var_temp_iter){
                            int col_ID = pair_temp.first;
                            int var_ID_temp = this->solver.locals[node_iter].ID.variable[var_temp_iter];

                            if(col_ID == var_ID_temp){
                                this->solver.locals[node_iter].infos.constraints[constr_iter].prime.neighbor.push_back(pair_temp);

                                var_temp_now = var_temp_iter + 1;

                                break_flag = 1;
                                break;
                            }
                        }

                        // No break occur in 2nd loop = outside
                        if(!break_flag){
                            this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside.push_back(pair_temp);
                        }
                    }
                }
            }
        }

        void print_node_infos(){
            std::cout << "Print infos of each node.\n\n";

            for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                std::cout << "=======================================================\n";
                std::cout << "Node " << node_iter << "\n";
                std::cout << "=======================================================\n";

                int num_var_temp = this->solver.locals[node_iter].ID.variable.size();
                int num_constr_temp = this->solver.locals[node_iter].ID.constraint.size();

                std::cout << "Variables:\n";
                std::cout << "ID:\t";
                for(int var_iter = 0; var_iter < num_var_temp; ++ var_iter){
                    std::cout << this->solver.locals[node_iter].ID.variable[var_iter] << "\t";
                }
                std::cout << "\n\n";

                std::cout << "Constraints:\n";
                std::cout << "ID:\t";
                for(int constr_iter = 0; constr_iter < num_constr_temp; ++ constr_iter){
                    std::cout << this->solver.locals[node_iter].ID.constraint[constr_iter] << "\t";
                }
                std::cout << "\n\n";

                for(int var_iter = 0; var_iter < num_var_temp; ++ var_iter){
                    std::cout << "------------------------------------------\n";
                    std::cout << "Var " << this->solver.locals[node_iter].ID.variable[var_iter] << "\n";
                    std::cout << "------------------------------------------\n";
                    std::cout << "Prime:\n";
                    std::cout << "Self:\t" << this->solver.locals[node_iter].infos.variables[var_iter].prime.self << "\n";

                    std::cout << "Neighbor:\t";
                    for(int prime_iter = 0; prime_iter < this->solver.locals[node_iter].infos.variables[var_iter].prime.neighbor.size(); ++ prime_iter){
                        std::cout << "(" << this->solver.locals[node_iter].infos.variables[var_iter].prime.neighbor[prime_iter].first << ", " << this->solver.locals[node_iter].infos.variables[var_iter].prime.neighbor[prime_iter].second << ")\t";
                    }
                    std::cout << "\n";

                    std::cout << "Outside:\t";
                    for(int prime_iter = 0; prime_iter < this->solver.locals[node_iter].infos.variables[var_iter].prime.outside.size(); ++ prime_iter){
                        std::cout << "(" << this->solver.locals[node_iter].infos.variables[var_iter].prime.outside[prime_iter].first << ", " << this->solver.locals[node_iter].infos.variables[var_iter].prime.outside[prime_iter].second << ")\t";
                    }
                    std::cout << "\n\n";

                    std::cout << "Dual:\n";
                    std::cout << "Neighbor:\t";
                    for(int dual_iter = 0; dual_iter < this->solver.locals[node_iter].infos.variables[var_iter].dual.neighbor.size(); ++ dual_iter){
                        std::cout << "(" << this->solver.locals[node_iter].infos.variables[var_iter].dual.neighbor[dual_iter].first << ", " << this->solver.locals[node_iter].infos.variables[var_iter].dual.neighbor[dual_iter].second << ")\t";
                    }
                    std::cout << "\n";
                    std::cout << "Outside:\t";
                    for(int dual_iter = 0; dual_iter < this->solver.locals[node_iter].infos.variables[var_iter].dual.outside.size(); ++ dual_iter){
                        std::cout << "(" << this->solver.locals[node_iter].infos.variables[var_iter].dual.outside[dual_iter].first << ", " << this->solver.locals[node_iter].infos.variables[var_iter].dual.outside[dual_iter].second << ")\t";
                    }
                    std::cout << "\n\n";
                }

                for(int constr_iter = 0; constr_iter < num_constr_temp; ++ constr_iter){
                    std::cout << "------------------------------------------\n";
                    std::cout << "Constr " << this->solver.locals[node_iter].ID.constraint[constr_iter] << "\n";
                    std::cout << "------------------------------------------\n";
                    std::cout << "Prime:\n";

                    std::cout << "Neighbor:\t";
                    for(int prime_iter = 0; prime_iter < this->solver.locals[node_iter].infos.constraints[constr_iter].prime.neighbor.size(); ++ prime_iter){
                        std::cout << "(" << this->solver.locals[node_iter].infos.constraints[constr_iter].prime.neighbor[prime_iter].first << ", " << this->solver.locals[node_iter].infos.constraints[constr_iter].prime.neighbor[prime_iter].second << ")\t";
                    }
                    std::cout << "\n";

                    std::cout << "Outside:\t";
                    for(int prime_iter = 0; prime_iter < this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside.size(); ++ prime_iter){
                        std::cout << "(" << this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside[prime_iter].first << ", " << this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside[prime_iter].second << ")\t";
                    }
                    std::cout << "\n\n";
                }
            }
        }

        // Price gap set (rho must be fixed!!)
        void price_gap_set(double rho){
            for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                for(int var_iter = 0; var_iter < this->solver.locals[node_iter].ID.variable.size(); ++ var_iter){
                    int var_ID = this->solver.locals[node_iter].ID.variable[var_iter];
                    double coeff = this->solver.locals[node_iter].infos.variables[var_iter].prime.self;

                    this->obj.cost_funcs[var_ID].moc.gap = this->obj.cost_funcs[var_ID].moc.price / rho;
                    this->obj.cost_funcs[var_ID].moc.gap += coeff * this->obj.cost_funcs[var_ID].moc.quantity;
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
                // Local optimizer
                for(int node_iter = 0; node_iter < this->statistic.num_node; ++ node_iter){
                    int num_var_temp = this->solver.locals[node_iter].ID.variable.size();
                    int num_constr_temp = this->solver.locals[node_iter].ID.constraint.size();

                    // Calculate the constant sums for prime variables
                    Eigen::VectorXd var_prime_outside_sum = Eigen::VectorXd::Zero(num_var_temp);
                    Eigen::VectorXd var_dual_outside_sum = Eigen::VectorXd::Zero(num_var_temp);
                    for(int var_iter = 0; var_iter < num_var_temp; ++ var_iter){
                        for(int prime_iter = 0; prime_iter < this->solver.locals[node_iter].infos.variables[var_iter].prime.outside.size(); ++ prime_iter){
                            int var_ID = this->solver.locals[node_iter].infos.variables[var_iter].prime.outside[prime_iter].first;
                            double coeff = this->solver.locals[node_iter].infos.variables[var_iter].prime.outside[prime_iter].second;
                            var_prime_outside_sum(var_iter) += coeff * this->solver.sol.prime.variables.prev(var_ID);
                        }

                        for(int dual_iter = 0; dual_iter < this->solver.locals[node_iter].infos.variables[var_iter].dual.outside.size(); ++ dual_iter){
                            int constr_ID = this->solver.locals[node_iter].infos.variables[var_iter].dual.outside[dual_iter].first;
                            double coeff = this->solver.locals[node_iter].infos.variables[var_iter].dual.outside[dual_iter].second;
                            var_dual_outside_sum(var_iter) += coeff * (this->solver.sol.dual.variables.prev(constr_ID) - this->solver.constraint.boundary(constr_ID));
                        }
                    }

                    // Calculate the constant sums for dual variables
                    Eigen::VectorXd constr_prime_outside_sum = Eigen::VectorXd::Zero(num_constr_temp);
                    for(int constr_iter = 0; constr_iter < num_constr_temp; ++ constr_iter){
                        for(int prime_iter = 0; prime_iter < this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside.size(); ++ prime_iter){
                            int var_ID = this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside[prime_iter].first;
                            double coeff = this->solver.locals[node_iter].infos.constraints[constr_iter].prime.outside[prime_iter].second;
                            constr_prime_outside_sum(constr_iter) += coeff * this->solver.sol.prime.variables.prev(var_ID);
                        }
                    }

                    while(true){
                        // Optimize prime variables
                        for(int var_iter = 0; var_iter < num_var_temp; ++ var_iter){
                            int var_ID = this->solver.locals[node_iter].ID.variable[var_iter];

                            // Find sum of inner product of neighbor and outside terms
                            double inner_prod = var_prime_outside_sum(var_iter);

                            double var_prime_neighbor_sum = 0.;
                            for(int prime_iter = 0; prime_iter < this->solver.locals[node_iter].infos.variables[var_iter].prime.neighbor.size(); ++ prime_iter){
                                auto pair_temp = this->solver.locals[node_iter].infos.variables[var_iter].prime.neighbor[prime_iter];
                                var_prime_neighbor_sum += pair_temp.second * this->solver.sol.prime.variables.now(pair_temp.first);
                            }

                            double var_dual_neighbor_sum = 0.;
                            for(int dual_iter = 0; dual_iter < this->solver.locals[node_iter].infos.variables[var_iter].dual.neighbor.size(); ++ dual_iter){
                                auto pair_temp = this->solver.locals[node_iter].infos.variables[var_iter].dual.neighbor[dual_iter];
                                var_dual_neighbor_sum += pair_temp.second * (this->solver.sol.dual.variables.now(pair_temp.first) - this->solver.constraint.boundary(pair_temp.first));
                            }

                            inner_prod += var_prime_neighbor_sum + var_dual_neighbor_sum;

                            // Bisection method for finding KKT point
                            Eigen::Vector2i gap_price_ID(0, this->obj.cost_funcs[var_ID].moc.price.size() - 1);
                            int mid_price_ID = gap_price_ID.sum() / 2;
                            while(gap_price_ID(1) - gap_price_ID(0) > 1){
                                Eigen::Vector2d gap_rent;
                                gap_rent(0) = this->obj.cost_funcs[var_ID].moc.gap(gap_price_ID(0)) + inner_prod;
                                gap_rent(1) = this->obj.cost_funcs[var_ID].moc.gap(gap_price_ID(1)) + inner_prod;

                                double mid_rent;
                                mid_rent = this->obj.cost_funcs[var_ID].moc.gap(mid_price_ID) + inner_prod;

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
                            double sensitivity = this->solver.locals[node_iter].infos.variables[var_iter].prime.self;
                            if(gap_price_ID(0) % 2 == 1){
                                // price fixed case
                                this->solver.sol.prime.price_margin(var_ID) = this->obj.cost_funcs[var_ID].moc.price(mid_price_ID);
                                this->solver.sol.prime.variables.now(var_ID) = -this->solver.sol.prime.price_margin(var_ID) / rho;
                                this->solver.sol.prime.variables.now(var_ID) -= inner_prod;
                                this->solver.sol.prime.variables.now(var_ID) /= sensitivity;
                            }
                            else{
                                // quantity fixed case
                                this->solver.sol.prime.variables.now(var_ID) = this->obj.cost_funcs[var_ID].moc.quantity(mid_price_ID);
                                this->solver.sol.prime.price_margin(var_ID) = -rho * (sensitivity * this->obj.cost_funcs[var_ID].moc.quantity(mid_price_ID) + inner_prod);
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
    };

    // Functions
    void radial_line_problem_set(opf_struct&, int, int, std::complex<double>, double, double, double);
}
