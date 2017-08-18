#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "Adaptive.hh"
#include "Compare.hh"

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
    double H[9]={ 1.6, 0, 0,
                  0, 0.8, 0,
                  0, 0, .1};
    double c[3] = {9.3,1,0};
    G->addEllipsoid(H, c, INNER);
};

auto unicycleAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
//    double h1[4] = {-1,1.2,-0, 9};
//    O->addPolytope(4,H,h1, OUTER);

    double h2[4] = {-2.2, 2.3, -0, 4.3};
    O->addPolytope(4, H, h2, OUTER);

    double h3[4] = {-2.2, 2.3, -6.2, 10};
    O->addPolytope(4, H, h3, OUTER);

    double h4[4] = {-3.4, 3.5, -0, 8.2};
    O->addPolytope(4, H, h4, OUTER);

    double h5[4] = {-4.6, 4.8, -1.5, 10};
    O->addPolytope(4,H,h5, OUTER);

//    double h6[4] = {-5.8,6,-0,6};
//    O->addPolytope(4,H,h6, OUTER);

//    double h7[4] = {-5.8,6,-7,10};
//    O->addPolytope(4,H,h7, OUTER);

//    double h8[4] = {-7,7.2,-1,10};
//    O->addPolytope(4,H,h8, OUTER);

    double h9[4] = {-7.4, 7.6, -0, 8.2};
    O->addPolytope(4,H,h9, OUTER);

//    double h10[4] = {-8.4,9.3,-8.3,8.5};
//    O->addPolytope(4,H,h10, OUTER);

//    double h11[4] = {-9.3,10,-7.1,7.3};
//    O->addPolytope(4,H,h11, OUTER);

//    double h12[4] = {-8.4,9.3,-5.9,6.1};
//    O->addPolytope(4,H,h12, OUTER);

//    double h13[4] = {-9.3,10 ,-4.7,4.9};
//    O->addPolytope(4,H,h13, OUTER);

//    double h14[4] = {-8.4,9.3,-3.5,3.7};
//    O->addPolytope(4,H,h14, OUTER);

//    double h15[4] = {-9.3,10 ,-2.3,2.5};
//    O->addPolytope(4,H,h15, OUTER);
};

auto unicycleAddI = [](SymbolicSet* I) -> void {
    double q[3] = {0.5, 0.5, M_PI/2};
    I->addPoint(q);
};

void sub(double* lbX, double* ubX, double* etaX, double tau, double* lbU, double* ubU, double* etaU,
         double* etaRatio, double tauRatio, int numAbs, int nint, int readAb) {
    Compare<X_type, U_type> comp(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readAb, "scots.txt");
    comp.initializeReach(unicycleAddG, unicycleAddI, unicycleAddO);
    comp.computeAbstractions(sysNext, radNext);
    int earlyBreak = 1;
    comp.reachSCOTS(earlyBreak);
}

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {10, 10, M_PI+0.4};

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    int nint = 5;

    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    double etaRatio[dimX] = {3, 3, 3};
    double tauRatio = 3;

    int numAbs = 2;
    int startAbs = 1;
    int readXX = 0;
    int readAbs = 0;

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readXX, readAbs, "adaptive.txt");
    abs.initializeReach(unicycleAddG, unicycleAddI, unicycleAddO);
    abs.computeAbstractions(sysNext, radNext);

    int minToGoCoarser = 1;
    int minToBeValid = 5;
    int earlyBreak = 1;

    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
    sub(lbX, ubX, etaX, tau, lbU, ubU, etaU, etaRatio, tauRatio, numAbs, nint, readAbs);
}



