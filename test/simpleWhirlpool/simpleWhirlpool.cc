#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "AlwaysEventually.hh"
#include "Compare.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define dimX 2
#define dimU 2

const double w[dimX] = {0.05, 0.05};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* we integrate the whirlpool ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](X_type &xx, const X_type &x, U_type &u) -> void {
        xx[0] = (-0.5 * x[0] + 1 * x[1]) * u[1];
        xx[1] = (-0.5 * x[0] + 0.5 * x[1] + u[0]) * u[1];
    };
    solver(ODE, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {        
        drdt[0] = -0.25*r[0] + 1*r[1] + w[0];
        drdt[1] = -0.25*r[0] + 0.5*r[1] + w[1];
    };
    solver(radODE, r, u);
};


auto whirlpoolAddO = [](SymbolicSet* O) -> void {
    ;
};

auto whirlpoolAddG = [](SymbolicSet* G) -> void {
    double H[dimX*2] = {1, 0,
                  0, 0.5};
    double c[dimX] = {0, 4};
    G->addEllipsoid(H, c, INNER);
};

auto whirlpoolAddI = [](SymbolicSet* I) -> void {
    double q[dimX] = {1, 1};
    I->addPoint(q);
};

void sub(double* lbX, double* ubX, double* etaX, double tau, double* lbU, double* ubU, double* etaU,
         double* etaRatio, double tauRatio, int numAbs, int nint, int readAb, int earlyBreak) {
    Compare<X_type, U_type> comp(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readAb, "scots.txt");
    comp.initializeReach(whirlpoolAddG, whirlpoolAddI, whirlpoolAddO);
    comp.computeAbstractions(sysNext, radNext);
//    comp.reachSCOTS(earlyBreak);
    comp.alwaysEventuallySCOTS();
}

int main() {

    double lbX[dimX]={-6, -6};
    double ubX[dimX]={ 6,  6};

    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2,   1};
    double etaU[dimU]= {0.5, 0.2};

    double etaRatio[dimX] = {3, 3};
    double tauRatio = 3;
    int nint = 5;

    double etaX[dimX]= {0.6, 0.6};
    double tau = 0.9;
    int numAbs = 2;

    int readXX = 0; // if X or U has changed, needs to be 0
    int readAbs = 0; // if X or U or O or dynamics has changed, needs to be 0


    AlwaysEventually<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readXX, readAbs, "adaptive.txt");

    abs.initialize(whirlpoolAddO);
    abs.initializeReach(whirlpoolAddG, whirlpoolAddI);
    abs.computeAbstractions(sysNext, radNext);

    int startAbs = 1;
    int minToGoCoarser = 2;
    int minToBeValid = 5;

//    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
    abs.alwaysEventually(startAbs, minToGoCoarser, minToBeValid);
//    sub(lbX, ubX, etaX, tau, lbU, ubU, etaU, etaRatio, tauRatio, numAbs, nint, readAbs, earlyBreak);

}
