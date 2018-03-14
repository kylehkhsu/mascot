#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

//#include "ReachAndStay.hh"
#include "AdaptAbsSafe.hh"
//#include "Compare.hh"

using namespace std;
using namespace scots;
using namespace helper;

/* state space dim */
#define dimX 2
#define dimU 1

/* disturbance */
const double w[1] = {0};

/* angular speed of the spiral (fixed) */
const double omega = 0.1;

/* maximum radiun */
const double maxRadius = 2;

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* system ODE (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto sysODE = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
        dxdt[0] = -x[0] + u[0]; /* radius in polar coordinate */
        dxdt[1] = omega; /* angular velocity */
    };
    solver(sysODE, x, u);
    if (x[0] > maxRadius) {
        x[0] = maxRadius;
    }
    if (x[0] < 0) {
        x[0] = -x[0];
        x[1] = x[1] + M_PI;
    }
    x[1] = fmod(x[1],(2*M_PI));
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

// Safe set
auto spiralAddS = [](SymbolicSet* S) -> void {
	/* avoid the obstacle */
    S->addGridPoints();
	double H[4 * 2] = { -1, 0,
						1, 0,
                        0, -1,
						0, 1};
	double c1[4] = { -1.5, 2, -M_PI + 0.3, M_PI + 0.3};
//    S->remPolytope(4, H, c1, INNER);
};


//void sub(double* lbX, double* ubX, double* etaX, double tau, double* lbU, double* ubU, double* etaU,
//         double* etaRatio, double tauRatio, int numAbs, int nint, int readAb) {
//    Compare<X_type, U_type> comp(dimX, lbX, ubX, etaX, tau,
//                                 dimU, lbU, ubU, etaU,
//                                 etaRatio, tauRatio, nint,
//                                 numAbs, readAb, "scots.txt");
//    comp.initializeReach(pendulumAddG, pendulumAddI, pendulumAddO);
//    comp.computeAbstractions(sysNext, radNext);
//    int earlyBreak = 1;
//    comp.reachSCOTS(earlyBreak);
//}

int main() {

    double lbX[dimX]={0, M_PI};
    double ubX[dimX]={maxRadius+0.1, 2*M_PI};

    double lbU[dimU]={-0.5};
    double ubU[dimU]= {0.5};
    double etaU[dimU]= {0.05};

    double etaRatio[dimX] = {2, 2};
    double tauRatio = 2;
    int nint = 5;

    double etaX[dimX]= {0.05, 0.1};
    double tau = 0.02;
    int numAbs = 1;

	int p = 2;
	
	X_type x;
	U_type u;

    //int readXX = 0; // if specification has changed, needs to be 0
    //int readAbs = 0;

    System spiral(dimX, lbX, ubX, etaX, tau,
                  dimU, lbU, ubU, etaU,
                  etaRatio, tauRatio, nint, numAbs);

	AdaptAbsSafe abs("spiral_3A.log");

    abs.initialize(&spiral, spiralAddS);

	TicToc timer;
	timer.tic();
	abs.onTheFlySafe(sysNext, radNext, x, u);
	clog << "-----------------------------------------------Total time: " << timer.toc() << " seconds.\n";
    /*abs.initializeSafe(pendulumAddG);
    abs.computeAbstractions(sysNext, radNext);*/

//    int startAbs = 1;
//    int minToGoCoarser = 2;
//    int minToBeValid = 3;
//    int earlyBreak = 1;

//    abs.reachAndStay(pendulumAddG, pendulumAddI, startAbs, minToGoCoarser, minToBeValid, earlyBreak);
    /*abs.safe();*/
}
