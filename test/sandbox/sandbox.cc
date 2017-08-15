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
using namespace helper;

/* data types for the ode solver */
typedef std::array<double,2> X_type;
typedef std::array<double,2> U_type;

auto simpleAddG = [](SymbolicSet* G) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {-7, 10, -7, 10};
    G->addPolytope(4, H, c, INNER);
};

int main() {
    int dimX = 2;
    int dimU = 2;

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

    abs.testProjections(simpleAddG);

}
