#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsReach.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define dimX 6
#define dimU 3

#define M3l -3
#define M3u 3
#define M4l 0
#define M4u 5
#define M5l 0
#define M5u 5

const double w[dimX] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.05};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* we integrate the ship ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, OdeSolver solver) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
    xx[0] = x[3]*cos(x[2])-x[4]*sin(x[2]);
    xx[1] = x[3]*sin(x[2]) + x[4]*cos(x[2]);
    xx[2] = x[5];
    xx[3] = 0.011*u[0] - 0.075*x[3];
    xx[4] = 0.01*u[1] - 0.001*u[2] - 0.38*x[4] - 0.98*x[5]*x[3] - 0.007*x[5];
    xx[5] = 0.045*u[2] - 0.001*u[1] - 0.077*x[4] - 0.013*x[5]*x[3] - 0.86*x[5];
  };
  solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */

auto radNext = [](X_type &r, U_type &u, OdeSolver solver) -> void {
  auto rhs=[](X_type& rr, const X_type& r, U_type& u) -> void {
    r[0] = ((M3u^2 + M4u^2)^0.5)*r[2] + r[3] + r[4];
    r[1] = ((M3u^2 + M4u^2)^0.5)*r[2] + r[3] + r[4];
    r[2] = r[5];
    r[3] = -0.075*r[3];
    r[4] = 0.98*M5u*r[3] - 0.38*r[4] + 0.007 + 0.98*M3u*r[5];
    r[5] = 0.013*M5u*r[3] + 0.077*r[4] - (0.86+0.013*M3l)*r[5];
  };
  solver(rhs,r,u);
};

auto shipAddG = [](SymbolicSet* G) -> void {
    double H[6*6]={-1, 0, 0, 0, 0, 0,
                    1, 0, 0, 0, 0, 0,
                    0,-1, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0,
                    0, 0, -1, 0, 0, 0,
                    0, 0, 1, 0, 0, 0
    };
  double h1[4] = {-7, 10, 0, 6.5, -M_PI/3, 2*M_PI/3};
    G->addPolytope(6, H, h1, INNER);
};

auto shipAddO = [](SymbolicSet* O) -> void {
    double H[6*6]={-1, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0,
                      0,-1, 0, 0, 0, 0,
                      0, 1, 0, 0, 0, 0,
                      0, 0, -1, 0, 0, 0,
                      0, 0, 1, 0, 0, 0
      };

    double h1[4] = {-2, 2.5, 0, 3, M_PI, M_PI};
    O->addPolytope(4, H, h1, OUTER);

  double h2[4] = {-5, 5.5, -3.5, 6.5, M_PI, M_PI};
    O->addPolytope(4, H, h2, OUTER);
};

int main() {

  double lbX[dimX] = {0, 0, -M_PI-0.4, M3l, M4l, M5l};
  double ubX[dimX] = {10, 6.5, M_PI+0.4, M3u, M4u, M5u};

  double lbU[dimU] = {0, -0.05, -0.1};
  double ubU[dimU] = {0.18, 0.05, 0.1};
  double etaU[dimU] = {.36, .2, 0.4};

    int nSubInt = 5;

    double etaX[dimX] = {1, 1, 1, 1, 1, 1};
    double tau = 1;

    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;
    int p = 2;
    int verbose = 1;

    int numAbs = 4;

    X_type x;
    U_type u;

    System ship(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    AdaptAbsReach abs("ship.log",verbose);
    abs.initialize(&ship, shipAddO, shipAddG);

    TicToc timer;
    timer.tic();
    abs.onTheFlyReach(p, sysNext, radNext, x, u);
    clog << "------------------------------------Total time: " << timer.toc() << " seconds.\n";
}
