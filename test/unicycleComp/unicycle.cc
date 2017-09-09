#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "Reach.hh"
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
    double H[9]={ 1.6,   0,   0,
                    0, 0.8,   0,
                    0,   0, 0.1};
    double c[3] = {8.1, 7.8, 0};
    G->addEllipsoid(H, c, INNER);

//    double H[9]={   1.6,  0,   0,
//                    0,  1.6,   0,
//                    0,  0, 0.1};
//    double c[3] = {4.5, 4.5, 0};
//    G->addEllipsoid(H, c, INNER);
};

auto unicycleAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double h1[4] = {-1.9, 2.3, -2, 10};
    O->addPolytope(4, H, h1, OUTER);

    double h2[4] = {-4.3, 4.7, -0, 7.1};
    O->addPolytope(4, H, h2, OUTER);

    double h3[4] = {-6.7, 7.1, -2, 10};
    O->addPolytope(4, H, h3, OUTER);

    double h4[4] = {-2.5, 3.2, -3.7, 4.6};
    O->addPolytope(4, H, h4, OUTER);

    double h5[4] = {-5.6, 6.5, -4.9, 6.5};
    O->addPolytope(4, H, h5, OUTER);

//    double h[4] = {-0.6, 5, -1.81, 2.39};
//    O->addPolytope(4, H, h, OUTER);
};

auto unicycleAddI = [](SymbolicSet* I) -> void {
    double q[3] = {1, 8.7, -M_PI/2};
    I->addPoint(q);

//    double q[3] = {4.5, 0.5, M_PI};
//    I->addPoint(q);
};

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {9.1, 9.1, M_PI+0.4};

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    int nSubInt = 5;

    double etaX[dimX] = {0.6/2/2, 0.6/2/2, 0.3/2/2};
    double tau = 0.9/2/2;

    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;

    int numAbs = 1;
    int startAbs = 0;
    int readXX = 0;
    int readAbs = 0;

    X_type x;
    U_type u;

    System unicycle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    Reach abs("unicycle1A.txt");
    abs.initialize(&unicycle, readXX, readAbs, unicycleAddO);
    abs.initializeReach(unicycleAddG, unicycleAddI);
    abs.computeAbstractions(sysNext, radNext, x, u);

    int minToGoCoarser = 5;
    int minToBeValid = 5;
    int earlyBreak = 1;

    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
}



