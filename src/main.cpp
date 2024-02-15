// Decentralized optimization using ADMM, relaxed constraints
// DC Power Flow Assumption, ignoring reactive power flow

#include "header.h"

int main()
{
    int num_sub = 2;
    int num_node = 6;
    int num_price = 10;
    std::complex<double> y_l (0., .5);
    double theta_limit = ADMM::pi() / 6.;
    double I_limit = 2.;
    double total_load = 2.;
    double penalty_price_voltage = 1000.;

    ADMM::opf_structs opfs;
    radial_line_problem_split_set(opfs, num_sub, num_node, num_price, y_l, theta_limit, I_limit, total_load, penalty_price_voltage);
//    opfs.solve_root_many_iteration(10, 1E-6, 1E-6, I_limit, 1);
    opfs.solve_root_many(1E-6, 1E-6, 0.6, 1);

//    ADMM::opf_struct opf;
//    radial_line_problem_set(opf, num_node, num_price, y_l, theta_limit, I_limit, total_load, std::numeric_limits<double>::infinity());
//    opf.solve_root(1E-6, 1E-6, 1);

    return 0;
}
