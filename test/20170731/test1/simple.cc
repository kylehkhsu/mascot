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
auto sysNext = [](X_type &x, U_type &u, double tau) -> void {
    x[0] += u[0] * tau;
    x[1] += u[1] * tau;
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau) -> void {
    r[0] = r[0] + abs(u[0]) * tau / 2;
    r[1] = r[1] + abs(u[1]) * tau / 2;
};

auto simpleAddO = [](SymbolicSet* O) -> void {

    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double h1[4] = {1, 14, -3, 5};
    O->addPolytope(4, H, h1, OUTER);
    double h2[4] = {-6, 21, -10, 12};
    O->addPolytope(4, H, h2, OUTER);
    double h3[4] = {-4, 16, -6, 9};
    O->addPolytope(4, H, h3, OUTER);
};

auto simpleAddG = [](SymbolicSet* G) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {-16, 19, -16, 19};
    G->addPolytope(4, H, c, INNER);
//    double c2[4] = {-12, 18, -17, 18};
//    G->addPolytope(4, H, c2, INNER);
};

auto simpleAddI = [](SymbolicSet* I) -> void {
    double q[2] = {1, 1};
    I->addPoint(q);
};

int main() {
    TicToc tt;
    Cudd mgr;
    double lbX[dimX]={0, 0};
    double ubX[dimX]={20, 20};
    double etaX[dimX]= {0.6, 0.6};
    double lbU[dimU]={-3, -3};
    double ubU[dimU]= {3, 3};
    double etaU[dimU]= {0.5, 0.5};
    double tau = 0.6;
    double etaRatio = 3;
    double tauRatio = 3;
    int numAbs = 2;
    int startAbs = 1;

    tt.tic();
    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio,
                                 numAbs, startAbs,
                                 simpleAddO, sysNext, radNext);

//    abs.reachBasicDebug(simpleAddG, simpleAddI, simpleAddO);
    abs.reach(simpleAddG, simpleAddI, simpleAddO);
    tt.toc();


//    abs.reachSCOTS();


}
