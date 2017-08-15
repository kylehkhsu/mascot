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
int dimX = 3+2;
int dimU = 2;

/* data types for the ode solver */
typedef std::array<double,3+2> X_type;
typedef std::array<double,2> U_type;

/* we integrate the unicycle ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
      xx[0] = u[0]*std::cos(x[2]);
      xx[1] = u[0]*std::sin(x[2]);
      xx[2] = u[1];
      xx[3] = 0;
      xx[4] = 0;
  };
  solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0]+r[2]*std::abs(u[0]) * tau;
    r[1] = r[1]+r[2]*std::abs(u[0]) * tau;
    r[2] = r[2];
    r[3] = r[3] / 1.001 * tau / 0.9;
    r[4] = r[4] / 1.001 * tau / 0.9;
};

auto unicycleAddG = [](SymbolicSet* G) -> void {
    auto f = [](double* x)->bool {
        return sqrt(pow(x[0]-x[3], 2) + pow(x[1]-x[4], 2)) <= 0.3;
    };
    G->addByFunction(f);
};

auto unicycleAddO = [](SymbolicSet* O) -> void {
    ;
};

auto unicycleAddI = [](SymbolicSet* I) -> void {
    double q[5] = {0.2, 0.2, -M_PI/2, 1.8, 1.8};
    I->addPoint(q);
};

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4, 0, 0};
    double ubX[dimX] = {2, 2, M_PI+0.4, 2, 2};

    double lbU[dimU] = {-1, -1.5};
    double ubU[dimU] = {1, 1.5};
    double etaU[dimU] = {1, .5};

    double etaX[dimX] = {0.6, 0.6, 0.3, 0.2, 0.2};
    double tau = 0.9;

    double etaRatio[dimX] = {3, 3, 3, 1, 1};
    double tauRatio = 3;
    int nint = 5;
    int numAbs = 2;
    int read = 0;


    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs,
                                 unicycleAddG, unicycleAddI, unicycleAddO,
                                 read);

    abs.computeAbstractions(sysNext, radNext, read);
	
    int startAbs = 1;
    int minToGoCoarser = 5;
    int minToBeValid = 5;
	
    abs.reach(startAbs, minToGoCoarser, minToBeValid, 1);
//    abs.reachSCOTS();

}



