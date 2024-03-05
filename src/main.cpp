// Decentralized optimization using ADMM, prime-dual gradient method
// DC Power Flow Assumption, ignoring reactive power flow

#include "header.h"

int main()
{
    int num_node = 111;
    int num_price = 10;
    std::complex<double> y_l (0., 1.);
    double theta_limit = ADMM::pi();
    double I_limit = 2.;
    double total_load = 2.;

    ADMM::opf_struct opf;
    radial_line_problem_set(opf, num_node, num_price, y_l, theta_limit, I_limit, total_load, std::numeric_limits<double>::infinity());
    opf.solve_ini();
    opf.solve_root(pow(10., -1.), pow(10., -1.), theta_limit, 1E7, 1, 1);
    opf.solve_root(pow(10., -2.), pow(10., -2.), theta_limit, 1E6, 1, 1);
    opf.solve_root(pow(10., -3.), pow(10., -3.), theta_limit, 1E5, 1, 1);
    opf.solve_root(pow(10., -4.), pow(10., -4.), theta_limit, 1E4, 1, 1);
    opf.solve_root(pow(10., -5.), pow(10., -5.), theta_limit, 1E3, 1, 1);
    opf.solve_root(pow(10., -6.), pow(10., -6.), theta_limit, 1E2, 1, 1);
    opf.solve_root(pow(10., -7.), pow(10., -7.), theta_limit, 1E1, 1, 1);
    opf.solve_root(pow(10., -9.), pow(10., -9.), theta_limit, 1E0, 1, 1);

    return 0;
}
