#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsReach.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define dimX 3
#define dimU 2

const double w[dimX] = {0.05, 0.05};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* we integrate the unicycle ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
      xx[0] = u[0]*std::cos(x[2]);
      xx[1] = u[0]*std::sin(x[2]);
      xx[2] = u[1];
  };
  solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */

auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0] + (r[2]*std::abs(u[0]) + w[0]) * tau;
    r[1] = r[1] + (r[2]*std::abs(u[0]) + w[1]) * tau;
};

auto unicycleAddG = [](SymbolicSet* G) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
    double h1[4] = {-10.79, 11.3, -0.1, 0.61};
    G->addPolytope(4, H, h1, INNER);
};

auto unicycleAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double h1[4] = {-1.9, 2.3, -1.9, 12};
    O->addPolytope(4, H, h1, OUTER);

    double h2[4] = {-4.3, 4.7, -0, 10.1};
    O->addPolytope(4, H, h2, OUTER);

    double h3[4] = {-6.7, 7.1, -1.9, 12};
    O->addPolytope(4, H, h3, OUTER);

    double h6[4] = {-9.1, 9.5, -0, 10.1};
    O->addPolytope(4, H, h6, OUTER);

    double h4[4] = {-2.5, 3.2, -3.7, 4.6};
    O->addPolytope(4, H, h4, OUTER);

    double h5[4] = {-5.39, 6.5, -4.9, 6.5};
    O->addPolytope(4, H, h5, OUTER);
    ;
};

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {11.4, 11.4, M_PI+0.4};

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    int nSubInt = 5;

    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;
    int p = 2;

    int numAbs = 3;

    X_type x;
    U_type u;

    System unicycle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    AdaptAbsReach abs("unicycle.log");
    abs.initialize(&unicycle, unicycleAddO, unicycleAddG);

    TicToc timer;
    timer.tic();
    abs.onTheFlyReach(p, sysNext, radNext, x, u);
    clog << "------------------------------------Total time: " << timer.toc() << " seconds.\n";
}



