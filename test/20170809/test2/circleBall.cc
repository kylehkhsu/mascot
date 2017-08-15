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
int dimX = 2+3;
int dimU = 2;

/* data types for the ode solver */
typedef std::array<double,2+3> X_type;
typedef std::array<double,2> U_type;

/* we integrate the circle ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](X_type &xx, const X_type &x, U_type &u) -> void {
        xx[0] = (-0.5 * x[0] + 1 * x[1]) * u[1];
        xx[1] = (-0.5 * x[0] + 0.5 * x[1] + u[0]) * u[1];
        xx[2] = 0;
        xx[3] = x[4];
        xx[4] = -(x[3]-3.5);
    };
    solver(ODE, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0];
    r[1] = r[1];
    r[2] = 0;
    r[3] = 0;
    r[4] = 0;
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
//    double H[4] = {1, 0,
//                   0, 1};
//    double c[2] = {0, 5};
//    G->addEllipsoid(H, c, INNER);

    auto f = [](double* x)->bool {
        return sqrt(pow(x[0]-x[2], 2) + pow(x[1]-x[3], 2)) <= 1;
    };

    G->addByFunction(f);

};

auto circleAddI = [](SymbolicSet* I) -> void {
    double q[5] = {1, 1, 0, 4.6, 0};
    I->addPoint(q);
};

int main() {

    double lbX[dimX]={-6, -6, -0.1, 2, -1.5};
    double ubX[dimX]={ 6,  6,  0.1, 5,  1.5};

    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2,   1};
    double etaU[dimU]= {0.5, 0.2};

    double etaRatio[dimX] = {3, 3, 1, 1, 1};
    double tauRatio = 3;
    int nint = 5;

    double etaX[dimX]= {0.6, 0.6, 0.2, 0.2, 0.2};
    double tau = 0.9;
    int numAbs = 2;

    int read = 1; // if X or U or O or dynamics has changed, needs to be 0


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
    int minToGoCoarser = 2;
    int minToBeValid = 5;
//    abs.reach(startAbs, minToGoCoarser, minToBeValid, 1);
    abs.alwaysEventually(startAbs, minToGoCoarser, minToBeValid, 1);




    //    abs.reachSCOTS(startAbs);
}
