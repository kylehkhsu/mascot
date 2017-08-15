#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "Adaptive.hh"

using namespace std;
using namespace scots;

/* state space dim */
int dimX = 2;
int dimU = 2;

/* disturbance */
const double w[2] = {0.3, 0.3};

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
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        double L[2][2];
        L[0][0] = 0;
        L[0][1] = 0;
        L[1][0] = 0;
        L[1][1] = 0;

        drdt[0] = L[0][0]*r[0] + L[0][1]*r[1] + w[0];
        drdt[1] = L[1][0]*r[0] + L[1][0]*r[1] + w[1];
    };
    solver(radODE, r, u);
};

auto simpleAddO = [](SymbolicSet* O) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double h1[4] = {-4, 5, 1, 9};
    O->addPolytope(4, H, h1, OUTER);
};

auto simpleAddG = [](SymbolicSet* G) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {-8, 10, -8, 10};
    G->addPolytope(4, H, c, INNER);
};

auto simpleAddI = [](SymbolicSet* I) -> void {
    double q[2] = {1, 1};
    I->addPoint(q);
};

int main() {

    double lbX[dimX]={0, 0};
    double ubX[dimX]={10, 10};

    double lbU[dimU]={-1, -1};
    double ubU[dimU]= {1, 1};
    double etaU[dimU]= {0.5, 0.5};

    double etaRatio[dimX] = {3, 3};
    double tauRatio = 3;
    int nint = 5;

    double etaX[dimX]= {0.6, 0.6};
    double tau = 0.9;
    int numAbs = 2;

    int readXX = 0; // if specification has changed, needs to be 0
    int readAbs = 0;

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readXX, readAbs);

    abs.initializeReach(simpleAddG, simpleAddI, simpleAddO);
    abs.computeAbstractions(sysNext, radNext);


    int startAbs = 1;
    int minToGoCoarser = 1;
    int minToBeValid = 2;
    abs.reach(startAbs, minToGoCoarser, minToBeValid, 1);

}
