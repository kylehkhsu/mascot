#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsReach.hh"
#include "Reach.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define dimX 3
#define dimU 2

const double w[dimX] = {0.05, 0.05};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {

  /* the ode describing the vehicle */
  auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
      double alpha=std::atan(std::tan(u[1])/2.0);
      xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
      xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
      xx[2] = u[0]*std::tan(u[1]);
  };
  solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */

auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    double c = std::abs(u[0]*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1));
    r[0] = r[0]+c*r[2]*tau;
    r[1] = r[1]+c*r[2]*tau;
};

auto vehicleAddG = [](SymbolicSet* G) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
    double h1[4] = {-2.5, 3.5, 1, 0.7};
    G->addPolytope(4, H, h1, INNER);
};

auto vehicleAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double h1[4] = {-2.1, 2.3, 0, 1.2};
    O->addPolytope(4, H, h1, OUTER);
    ;
};

auto vehicleAddI = [](SymbolicSet* I) -> void {
    double q[3] = {0.25, 0.25, 0};
    I->addPoint(q);
};

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {6.0, 1.8, M_PI+0.4};

    double lbU[dimU] = {-1, -1};
    double ubU[dimU] = {1, 1};
    double etaU[dimU] = {.3, .3};

    int nSubInt = 5;

    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;

    int numAbs = 3;
//    int startAbs = 0;
//    int readAbs = 0;

    X_type x;
    U_type u;

    System vehicle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    AdaptAbsReach abs("vehicle3A.txt");
    abs.initialize(&vehicle, vehicleAddO, vehicleAddG);

    TicToc tt_tot;
    tt_tot.tic();
    abs.onTheFlyReach(sysNext, radNext, x, u);
    clog << "------------------------------------Total time:";
    tt_tot.toc();

//    abs.initializeReach(vehicleAddG, vehicleAddI);
//    abs.computeAbstractions(sysNext, radNext, x, u);
//
//    int minToGoCoarser = 6;
//    int minToBeValid = 6;
//    int earlyBreak = 0;
//
//    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
}



