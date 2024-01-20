// Decentralized optimization using ADMM, relaxed constraints
// DC Power Flow Assumption, ignoring reactive power flow

#include "header.h"

int main()
{
    int num_node = 6;
    int num_price = 100;
    std::complex<double> y_l (0., 1.);
    double theta_limit = ADMM::pi() / 6.;
    double I_limit = 1.;
    double total_load = 2.;

    ADMM::opf_structs opfs;
    radial_line_problem_split_set(opfs, num_node, num_price, y_l, theta_limit, I_limit, total_load);
    opfs.solve_root_many(1E-6, 1E-6, 1);

//    ADMM::opf_struct opf;
//    radial_line_problem_set(opf, num_node, num_price, y_l, theta_limit, I_limit, total_load);
//    opf.solve_root(1E-6, 1E-6, 1);

    // Profile of algorith
//    for(int test_iter = 0; test_iter < 1; ++ test_iter){
//        if(test_iter % 10 == 0){
//            std::cout << test_iter << "\n";
//        }
//        opf.solve_RK4(1E-6, 1E-4, 1);
//        opf.solve_root(1E-6, 1E-4, 1);
//    }

    return 0;
}
