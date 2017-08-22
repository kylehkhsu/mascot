#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "ReachAndStay.hh"
#include "Compare.hh"

using namespace std;
using namespace scots;

/* state space dim */
#define dimX 4
#define dimU 1

/* disturbance */
const double w[1] = {0};

/* system parameters in S.I. units */
const double M = 2.4; // Mass of cart
const double m = 0.23; // Mass of bob
const double l = 0.36; // length of rod
const double fr = 0.1; // friction of the pendulum
const double g = 9.81; // acceleration due to gravity

/* data types for the ode solver */
typedef std::array<double,4> X_type;
typedef std::array<double,1> U_type;

/* system ODE (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto sysODE = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
        dxdt[0] = x[1];
        dxdt[1] = (-m*g*sin(x[2])*cos(x[2]) + m*l*pow(x[3],2)*sin(x[2]) + fr*m*x[3]*cos(x[2]) + u[0])/(M + (1-pow(cos(x[2]),2))*m);
        dxdt[2] = x[3];
        dxdt[3] = ((M+m)*(g*sin(x[2]) - fr*x[3]) - (l*m*pow(x[3],2)*sin(x[2])+u[0])*cos(x[2]))/(l*(M+ (1-pow(cos(x[2]),2))*m));
    };
    solver(sysODE, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
//    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
//        double L[4][4];
//        L[0][0] = 0;
//        L[0][1] = 1;
//        L[0][2] = 0;
//        L[0][3] = 0;
//        L[1][0] = 0;
//        L[1][1] = 0;
//        L[1][2] = abs(0.0399*u[0] + 1.02);
//        L[1][3] = 0.0882;
//        L[2][0] = 0;
//        L[2][1] = 0;
//        L[2][2] = 0;
//        L[2][3] = 1;
//        L[3][0] = 0;
//        L[3][1] = 0;
//        L[3][2] = abs(1.27*u[0] + 32.9);
//        L[3][3] = -0.209;
//
//        drdt[0] = L[0][0]*r[0] + L[0][1]*r[1] + L[0][2]*r[2] + L[0][3]*r[3];
//        drdt[1] = L[1][0]*r[0] + L[1][1]*r[1] + L[1][2]*r[2] + L[1][3]*r[3];
//        drdt[2] = L[2][0]*r[0] + L[2][1]*r[1] + L[2][2]*r[2] + L[2][3]*r[3];
//        drdt[3] = L[3][0]*r[0] + L[3][1]*r[1] + L[3][2]*r[2] + L[3][3]*r[3] + w[0];
//    };
//    solver(radODE, r, u);
    ;
};

auto pendulumAddO = [](SymbolicSet* O) -> void {
    ;
};

auto pendulumAddG = [](SymbolicSet* G) -> void {
    /* balance the pendulum within -0.09 rad = 2*pi-0.09 rad to 0.09 rad */
    double H[2*2]={-1, 0,
                    1, 0};
    double c1[2] = {0, 0.18};
    G->addPolytope(2, H, c1, OUTER);
    
    double c2[2] = {-2*M_PI+0.18, 2*M_PI};
    G->addPolytope(2, H, c2, OUTER);
};


//auto pendulumAddG = [](SymbolicSet* G) -> void {
//    /* balance the pendulum within -0.02 rad = 2*pi-0.02 rad to 0.02 rad near the center of the rail (+/- 0.1) */
//    double H[4*4]={-1, 0, 0, 0,
//                      1, 0, 0, 0,
//                      0, 0, -1, 0,
//                      0, 0, 1, 0};
//    double c1[4] = {0.1, 0.1, 0, 0.02};
//    G->addPolytope(4, H, c1, INNER);
//
//    double c2[4] = {0.1, 0.1, -2*M_PI+0.02, 2*M_PI};
//    G->addPolytope(4, H, c2, INNER);
//};

auto pendulumAddI = [](SymbolicSet* I) -> void {
    double q[4] = {0, 0, M_PI, 0};
    I->addPoint(q);
};

void sub(double* lbX, double* ubX, double* etaX, double tau, double* lbU, double* ubU, double* etaU,
         double* etaRatio, double tauRatio, int numAbs, int nint, int readAb) {
    Compare<X_type, U_type> comp(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readAb, "scots.txt");
    comp.initializeReach(pendulumAddG, pendulumAddI, pendulumAddO);
    comp.computeAbstractions(sysNext, radNext);
    int earlyBreak = 1;
    comp.reachSCOTS(earlyBreak);
}

int main() {

    double lbX[dimX]={-0.5, -5, 0, -1};
    double ubX[dimX]={0.5, 5, 2*M_PI, 1};

    double lbU[dimU]={-24};
    double ubU[dimU]= {24};
    double etaU[dimU]= {0.1};

    double etaRatio[dimX] = {4, 4, 2, 6};
    double tauRatio = 3;
    int nint = 5;

    double etaX[dimX]= {0.2, 2, 0.01, 0.6};
    double tau = 0.9;
    int numAbs = 2;

    int readXX = 0; // if specification has changed, needs to be 0
    int readAbs = 0;

    System system(dimX, lbX, ubX, etaX, tau,
                  dimU, lbU, ubU, etaU,
                  etaRatio, tauRatio, nint, numAbs);

    //ReachAndStay<X_type, U_type> abs("adaptive.txt");
    Safe<X_type, U_type> abs("adaptive.txt");

    abs.initialize(&system, readXX, readAbs, pendulumAddO);
    abs.initializeSafe(pendulumAddG);

    abs.computeAbstractions(sysNext, radNext);

//    int startAbs = 1;
//    int minToGoCoarser = 2;
//    int minToBeValid = 3;
//    int earlyBreak = 1;

//    abs.reachAndStay(pendulumAddG, pendulumAddI, startAbs, minToGoCoarser, minToBeValid, earlyBreak);
    abs.safe();
}
