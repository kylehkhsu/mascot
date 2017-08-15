#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "Adaptive.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

using namespace std;
using namespace scots;

/* state space dim */
int dimX = 2;
int dimU = 2;

/* data types for the ode solver */
typedef std::array<double,2> X_type;
typedef std::array<double,2> U_type;

/* we integrate the circle ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](X_type &xx, const X_type &x, U_type &u) -> void {
        xx[0] = (-1 * x[0] + 2 * x[1]) * u[1];
        xx[1] = (-1 * x[0] + 1 * x[1] + u[0]) * u[1];
    };
    solver(ODE, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0];
    r[1] = r[1];
};

auto circleAddO = [](SymbolicSet* O) -> void {
//    double H[4*2]={-1, 0,
//                      1, 0,
//                      0,-1,
//                      0, 1};
//    double h1[4] = {1, 3, -16, 18};
//    O->addPolytope(4, H, h1, OUTER);
//    double h2[4] = {-3.8, 5, -16, 20};
//    O->addPolytope(4, H, h2, OUTER);
//    double h3[4] = {-7, 13, -1.5, 20};
//    O->addPolytope(4, H, h3, OUTER);
};

auto circleAddG = [](SymbolicSet* G) -> void {
    double H[4] = {1, 0,
                   0, 1};
    double c[2] = {0, 4};
    G->addEllipsoid(H, c, INNER);
};

auto circleAddI = [](SymbolicSet* I) -> void {
    double q[2] = {1, 1};
    I->addPoint(q);
};

int main() {

    double lbX[dimX]={-7, -7};
    double ubX[dimX]={ 7,  7};

    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2, 1};
    double etaU[dimU]= {0.2, 0.2};

    double etaRatio[dimX] = {3, 3};
    double tauRatio = 3;
    int nint = 5;

    double etaX[dimX]= {0.6, 0.6};
    double tau = 0.9;
    int numAbs = 2;

    int read = 0; // if specification has changed, needs to be 0


    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs,
                                 circleAddG, circleAddI, circleAddO,
                                 read);
    //    abs.test(circleAddG);

    //    abs.wasteMemory();
    abs.computeAbstractions(sysNext, radNext, read);

    int startAbs = 1;
    int minToGoCoarser = 4;
    int minToBeValid = 4;
//    abs.reach(startAbs, minToGoCoarser, minToBeValid, 1);
    abs.alwaysEventually(startAbs, minToGoCoarser, minToBeValid, 1);




    //    abs.reachSCOTS(startAbs);
}
