#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "Reach.hh"
#include "Compare.hh"

using namespace std;
using namespace scots;

/* state space dim */
#define dimX 2
#define dimU 2

/* disturbance */
const double w[dimX] = {0.3, 0.3};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

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
    double c[4] = {-7.5, 10, -7.5, 10};
    G->addPolytope(4, H, c, INNER);
};

auto simpleAddI = [](SymbolicSet* I) -> void {
    double q[2] = {1, 1};
    I->addPoint(q);
};

void sub(double* lbX, double* ubX, double* etaX, double tau, double* lbU, double* ubU, double* etaU,
         double* etaRatio, double tauRatio, int numAbs, int nint, int readAb) {
    Compare<X_type, U_type> comp(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readAb, "scots.txt");
    comp.initializeReach(simpleAddG, simpleAddI, simpleAddO);
    comp.computeAbstractions(sysNext, radNext);
    int earlyBreak = 1;
    comp.reachSCOTS(earlyBreak);
}

//void testSynthesis() {
//    double lbX[dimX]={0, 0};
//    double ubX[dimX]={10, 10};

//    double lbU[dimU]={-1.3, -1.3};
//    double ubU[dimU]= {1.3, 1.3};
//    double etaU[dimU]= {0.5, 0.5};

//    double etaRatio[dimX] = {2, 2};
//    double tauRatio = 2;
//    int nint = 5;

//    double etaX[dimX]= {0.8, 0.8};
//    double tau = 1.2;
//    int numAbs = 3;

//    int readXX = 0; // if specification has changed, needs to be 0
//    int readAbs = 0;

//    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
//                                 dimU, lbU, ubU, etaU,
//                                 etaRatio, tauRatio, nint,
//                                 numAbs, readXX, readAbs, "adaptive.txt");
////    abs.testProjections(simpleAddG, 1);

//    abs.initializeReach(simpleAddG, simpleAddI, simpleAddO);
//    abs.computeAbstractions(sysNext, radNext);


//    int startAbs = 2;
//    int minToGoCoarser = 1;
//    int minToBeValid = 2;
//    int earlyBreak = 1;

//    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
//    sub(lbX, ubX, etaX, tau, lbU, ubU, etaU, etaRatio, tauRatio, numAbs, nint, readAbs);

//}

int main() {

    double lbX[dimX]={0, 0};
    double ubX[dimX]={10, 10};

    double lbU[dimU]={-1.3, -1.3};
    double ubU[dimU]= {1.3, 1.3};
    double etaU[dimU]= {0.5, 0.5};

    double etaRatio[dimX] = {2, 2};
    double tauRatio = 2;
    int nint = 5;

    double etaX[dimX]= {0.8, 0.8};
    double tau = 1.2;
    int numAbs = 3;

    int readXX = 0; // if specification has changed, needs to be 0
    int readAbs = 1;

    Reach<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                              dimU, lbU, ubU, etaU,
                              etaRatio, tauRatio, nint,
                              numAbs, readXX, readAbs, "adaptive.txt");
    abs.initialize(simpleAddO);
    abs.initializeReach(simpleAddG, simpleAddI);
    abs.computeAbstractions(sysNext, radNext);

    int startAbs = 2;
    int minToGoCoarser = 1;
    int minToBeValid = 2;
    int earlyBreak = 1;

    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);



}
