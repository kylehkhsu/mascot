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

#define dimX 2
#define dimU 2

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

int main() {
    Cudd mgr;

    double lbX[dimX] = {-3, -3};
    double ubX[dimX] = { 2,  2};
    double etaX[dimX] = {1, 1};
    double tau = 1;

    SymbolicSet X(mgr, dimX, lbX, ubX, etaX, tau);

    X.printInfo(1);

    const int minterm[6] = {0, 0, 0, 0, 0, 0};
    double p[2] = {0, 0};

    X.mintermToElement(minterm, p);
    printArray(p, 2);

    X.addPoint(p);
    X.printInfo(2);
    X.addGridPoints();

    X.writeToFile("X.bdd");

    SymbolicSet R(X);

    double H[4*2] = {-1,  0,
                      1,  0,
                      0, -1,
                      0,  1 };
    double h[4] = {2.5, 0.5, 2.5, 0.5};
    R.addPolytope(4, H, h, INNER);
//    R.addPolytope(4, H, h, OUTER);

    R.writeToFile("R.bdd");

    SymbolicSet E(X);
    double C[2*2] = {1, 0,
                     0, 1};
    double c[2] = {0, 0.5};

    E.addEllipsoid(C, c, OUTER);
    E.writeToFile("E.bdd");


}

