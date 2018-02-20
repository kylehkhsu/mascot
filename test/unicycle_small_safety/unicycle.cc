#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsSafe.hh"
//#include "Reach.hh"
//#include "UpfrontReach.hh"

using namespace std;
using namespace scots;

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

auto unicycleAddS = [](SymbolicSet* S) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double h1[4] = {-2.1, 2.3, 0, 1.2};
    S->addGridPoints();
    S->remPolytope(4, H, h1, OUTER);
    ;
};

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {6.0, 1.8, M_PI+0.4};

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    int nSubInt = 5;

    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;

    int numAbs = 3;
    int p = 2;

    X_type x;
    U_type u;

    System unicycle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    AdaptAbsSafe abs("unicycle_small_3A_safe.txt");
    abs.initialize(&unicycle, unicycleAddS);

    TicToc tt_tot;
    tt_tot.tic();
    abs.onTheFlySafe(sysNext, radNext, x, u);
    clog << "------------------------------------Total time:";
    tt_tot.toc();

//    int m = 2;
//    UpfrontReach abs("unicycle_small_3A_HSCC_recursive.log");
//    abs.initialize(&unicycle, 0, unicycleAddO);
//    abs.initializeReach(unicycleAddG, unicycleAddI);

//    TicToc timer;
//    timer.tic();
//    abs.computeAbstractions(sysNext, radNext, x, u);
//    abs.upfrontReach(m);
//    clog << "------------------------------------Total time: " << timer.toc() << " seconds.\n";
}



