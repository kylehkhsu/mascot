#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "Safe.hh"

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
        drdt[0] = (-0.5*r[0] + r[1]) * u[1] + w[0];
        drdt[1] = ( 0.5*r[0] + 0.5*r[1]) * u[1] + w[1];
    };
    solver(radODE, r, u);
};


auto whirlpoolAddO = [](SymbolicSet* O) -> void {
    double H[4]={ 0.8,   0,
                    0, 0.8};
    double c[2] = {4, 0};
    O->addEllipsoid(H, c, OUTER);
};

auto whirlpoolAddS = [](SymbolicSet* S) -> void {
    S->addGridPoints();
};

int main() {

    double lbX[dimX]={-6, -6};
    double ubX[dimX]={ 6,  6};

    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2,   1};
    double etaU[dimU]= {0.5, 0.2};

    double etaRatio[dimX] = {2, 2};
    double tauRatio = 2;
    int nint = 5;

    double etaX[dimX]= {0.6, 0.6};
    double tau = 0.6;
    int numAbs = 3;

    X_type x;
    U_type u;

    int readAbs = 0; // if X or U or O or dynamics has changed, needs to be 0

    System whirlpool(dimX, lbX, ubX, etaX, tau,
                dimU, lbU, ubU, etaU,
                etaRatio, tauRatio, nint, numAbs);

    Safe abs("whirlpoolSafe3Asafe2.txt");
    abs.initialize(&whirlpool, readAbs, whirlpoolAddO);
    abs.initializeSafe(whirlpoolAddS);
    abs.computeAbstractions(sysNext, radNext, x, u);
//    abs.safe();

    int itersToNextAbs = 100;
    abs.safe2(itersToNextAbs);

}
