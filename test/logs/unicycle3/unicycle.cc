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
    double c[3] = {9.3, 9, 0};
    G->addEllipsoid(H, c, INNER);
};

auto unicycleAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double h1[4] = {-1.9, 2.3, -2, 10};
    O->addPolytope(4, H, h1, OUTER);

    double h2[4] = {-4.3, 4.7, -0, 8.2};
    O->addPolytope(4, H, h2, OUTER);

    double h3[4] = {-6.7, 7.1, -2, 10};
    O->addPolytope(4, H, h3, OUTER);

    double E[9] = { 4, 0, 0,
                    0, 4, 0,
                    0, 0, 0};
    double c1[3] = {0.9, 6.9, 0};
    O->addEllipsoid(E, c1, OUTER);

    double c2[3] = {1, 3, 0};
    O->addEllipsoid(E, c2, OUTER);

    double c3[3] = {3.4, 2, 0};
    O->addEllipsoid(E, c3, OUTER);

    double c4[3] = {3, 7, 0};
    O->addEllipsoid(E, c4, OUTER);

    double c5[3] = {5.7, 6, 0};
    O->addEllipsoid(E, c5, OUTER);

    double c6[3] = {5.9, 2.5, 0};
    O->addEllipsoid(E, c6, OUTER);

    double c7[3] = {9.4, 1.8, 0};
    O->addEllipsoid(E, c7, OUTER);

    double c8[3] = {8.2, 6.2, 0};
    O->addEllipsoid(E, c8, OUTER);
};

auto unicycleAddI = [](SymbolicSet* I) -> void {
    double q[3] = {1, 9, -M_PI/2};
    I->addPoint(q);
};

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {10, 10, M_PI+0.4};

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    int nSubInt = 5;

    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    double etaRatio[dimX] = {3, 3, 3};
    double tauRatio = 3;

    int numAbs = 2;
    int startAbs = 1;
    int readXX = 1;
    int readAbs = 1;

    X_type x;
    U_type u;

    System unicycle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    Reach abs("unicycle2AinnerCoarserSymbolic.txt");
    abs.initialize(&unicycle, readXX, readAbs, unicycleAddO);
    abs.initializeReach(unicycleAddG, unicycleAddI);
    abs.computeAbstractions(sysNext, radNext, x, u);

    int minToGoCoarser = 1;
    int minToBeValid = 5;
    int earlyBreak = 1;

    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
}



