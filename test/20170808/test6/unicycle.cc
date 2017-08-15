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
int dimX = 3;
int dimU = 2;

/* data types for the ode solver */
typedef std::array<double,3> X_type;
typedef std::array<double,2> U_type;

/* we integrate the unicycle ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
      xx[0] = u[0]*std::cos(x[2]);
      xx[1] = u[0]*std::sin(x[2]);
      xx[2] = u[1];
  };
  solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0]+r[2]*std::abs(u[0]) * tau;
    r[1] = r[1]+r[2]*std::abs(u[0]) * tau;
};

auto unicycleAddG = [](SymbolicSet* G) -> void {
    double H[9]={ 1, 0, 0,
                  0, 0.5, 0,
                  0, 0, .1};
    double c[3] = {9,1,0};
    G->addEllipsoid(H, c, INNER);
};

auto unicycleAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
//    /* remove outer approximation of P={ x | H x<= h1 } form state space */
//    double h1[4] = {-1,1.2,-0, 9};
//    O->addPolytope(4,H,h1, OUTER);
//    /* remove outer approximation of P={ x | H x<= h2 } form state space */
//    double h2[4] = {-2.2,2.4,-0,5};
//    O->addPolytope(4,H,h2, OUTER);
//    /* remove outer approximation of P={ x | H x<= h3 } form state space */
//    double h3[4] = {-2.2,2.4,-6,10};
//    O->addPolytope(4,H,h3, OUTER);
//    /* remove outer approximation of P={ x | H x<= h4 } form state space */
//    double h4[4] = {-3.4,3.6,-0,9};
//    O->addPolytope(4,H,h4, OUTER);
//    /* remove outer approximation of P={ x | H x<= h5 } form state space */
    double h5[4] = {-4.6 ,4.8,-1,10};
    O->addPolytope(4,H,h5, OUTER);
//    /* remove outer approximation of P={ x | H x<= h6 } form state space */
    double h6[4] = {-5.8,6,-0,6};
    O->addPolytope(4,H,h6, OUTER);
//    /* remove outer approximation of P={ x | H x<= h7 } form state space */
    double h7[4] = {-5.8,6,-7,10};
    O->addPolytope(4,H,h7, OUTER);
//    /* remove outer approximation of P={ x | H x<= h8 } form state space */
//    double h8[4] = {-7,7.2,-1,10};
//    O->addPolytope(4,H,h8, OUTER);
//    /* remove outer approximation of P={ x | H x<= h9 } form state space */
//    double h9[4] = {-8.2,8.4,-0,8.5};
//    O->addPolytope(4,H,h9, OUTER);
//    /* remove outer approximation of P={ x | H x<= h10 } form state space */
//    double h10[4] = {-8.4,9.3,-8.3,8.5};
//    O->addPolytope(4,H,h10, OUTER);
//    /* remove outer approximation of P={ x | H x<= h11 } form state space */
//    double h11[4] = {-9.3,10,-7.1,7.3};
//    O->addPolytope(4,H,h11, OUTER);
//    /* remove outer approximation of P={ x | H x<= h12 } form state space */
//    double h12[4] = {-8.4,9.3,-5.9,6.1};
//    O->addPolytope(4,H,h12, OUTER);
//    /* remove outer approximation of P={ x | H x<= h13 } form state space */
//    double h13[4] = {-9.3,10 ,-4.7,4.9};
//    O->addPolytope(4,H,h13, OUTER);
//    /* remove outer approximation of P={ x | H x<= h14 } form state space */
//    double h14[4] = {-8.4,9.3,-3.5,3.7};
//    O->addPolytope(4,H,h14, OUTER);
//    /* remove outer approximation of P={ x | H x<= h15 } form state space */
//    double h15[4] = {-9.3,10 ,-2.3,2.5};
//    O->addPolytope(4,H,h15, OUTER);
};

auto unicycleAddI = [](SymbolicSet* I) -> void {
    double q[3] = {1, 8, -M_PI/2};
    I->addPoint(q);
};

int main() {
    TicToc tt;

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {10, 10, M_PI+0.4};

    double lbU[dimU] = {-1, -1.5};
    double ubU[dimU] = {1, 1.5};
    double etaU[dimU] = {.3, .2};

    int nint = 5;

    double etaX[dimX] = {0.6/3, 0.6/3, 0.3/3};
    double tau = 0.9/3;

    double etaRatio[dimX] = {3, 3, 3};
    double tauRatio = 3;

    int numAbs = 1;
    int startAbs = 0;
    int read = 0;


    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs,
                                 unicycleAddG, unicycleAddI, unicycleAddO,
                                 read);

    abs.computeAbstractions(sysNext, radNext, read);

	int minToGoCoarser = 5;
	int minToBeValid = 5;
	
//    abs.reach(startAbs, minToGoCoarser, minToBeValid);
    abs.reachSCOTS(startAbs);


}



