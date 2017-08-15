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
int dimX = 4;
int dimU = 2;

/* data types for the ode solver */
typedef std::array<double,4> X_type;
typedef std::array<double,2> U_type;

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    x[0] += u[0] * tau;
    x[1] += u[1] * tau;
    x[2] = x[2];
    x[3] = x[3];
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0] + abs(u[0]) * tau / 4;
    r[1] = r[1] + abs(u[1]) * tau / 4;
    r[2] = r[2] / 1.001 * tau / 2.7;
    r[3] = r[3] / 1.001 * tau / 2.7;
};

auto simpleAddO = [](SymbolicSet* O) -> void {
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
    ;
};

auto simpleAddG = [](SymbolicSet* G) -> void {

    auto f = [](double* x)->bool {
        return sqrt(pow(x[0]-x[2], 2) + pow(x[1]-x[3], 2)) <= 1;
    };

    G->addByFunction(f);
};

auto simpleAddI = [](SymbolicSet* I) -> void {
    double q[4] = {0.5, 0.5, 14, 14};
    I->addPoint(q);
};

int main() {

    double lbX[dimX]={0, 0, 0, 0};
    double ubX[dimX]={15, 15, 15, 15};
    double etaX[dimX]= {1.8, 1.8, 0.6, 0.6};
    double lbU[dimU]={-1, -1};
    double ubU[dimU]= {1, 1};
    double etaU[dimU]= {0.5, 0.5};
    double tau = 2.7;
    double etaRatio[dimX] = {3, 3, 1, 1};
    double tauRatio = 3;
    int nint = 5;
    int numAbs = 2;
    int read = 0; // if Xs or U or O or dynamics have changed, needs to be 0

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs,
                                 simpleAddG, simpleAddI, simpleAddO,
                                 read);
//    abs.test(simpleAddG);

//    abs.wasteMemory();
    abs.computeAbstractions(sysNext, radNext, read);

    int startAbs = 1;
    int minToGoCoarser = 3;
    int minToBeValid = 3;
    abs.reach(startAbs, minToGoCoarser, minToBeValid);

//    abs.reachSCOTS();
}
