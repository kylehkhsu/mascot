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

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    x[0] += u[0] * tau;
    x[1] += u[1] * tau;
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0] + abs(u[0]) * tau / 4;
    r[1] = r[1] + abs(u[1]) * tau / 4;
};

auto simpleAddO = [](SymbolicSet* O) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double h1[4] = {-8, 10, -8, 20};
    O->addPolytope(4, H, h1, OUTER);
    double h2[4] = {-10, 17, -8, 10};
    O->addPolytope(4, H, h2, OUTER);
    double h3[4] = {-12, 13, -14, 20};
    O->addPolytope(4, H, h3, OUTER);
    double h4[4] = {-14, 20, -13, 14};
    O->addPolytope(4, H, h4, OUTER);
};

auto simpleAddG = [](SymbolicSet* G) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {-15, 20, -15, 20};
    G->addPolytope(4, H, c, INNER);
};

auto simpleAddI = [](SymbolicSet* I) -> void {
    double q[2] = {1, 19};
    I->addPoint(q);
};

int main() {

    double lbX[dimX]={0, 0};
    double ubX[dimX]={20, 20};
    double etaX[dimX]= {1.8, 1.8};
    double lbU[dimU]={-3, -3};
    double ubU[dimU]= {3, 3};
    double etaU[dimU]= {0.5, 0.5};
    double tau = 1.8;
    double etaRatio[dimX] = {3, 3};
    double tauRatio = 3;
    int nint = 5;
    int numAbs = 3;
    int mode = 0; // if specification has changed, needs to be 0

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs,
                                 simpleAddG, simpleAddI, simpleAddO,
                                 mode);

//    abs.wasteMemory();
    abs.computeAbstractions(sysNext, radNext, mode);

    int startAbs = 2;
    int minToGoCoarser = 3;
    int minToBeValid = 3;
    abs.reach(startAbs, minToGoCoarser, minToBeValid);

//    abs.reachSCOTS();
}
