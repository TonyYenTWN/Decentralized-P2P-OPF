// Decentralized optimization using ADMM, relaxed constraints
// DC Power Flow Assumption, ignoring reactive power flow

#include "header.h"

int main()
{
    int num_node = 80;
    int num_price = 10;
    std::complex<double> y_l (0., 1.);
    double theta_limit = ADMM::pi() / 6.;
    double I_limit = 2.;
    double total_load = 2.;

    ADMM::opf_struct opf;
    radial_line_problem_set(opf, num_node, num_price, y_l, theta_limit, I_limit, total_load, std::numeric_limits<double>::infinity());
    opf.solve_root(1E-6, 1E-6, 1);

    return 0;
}
