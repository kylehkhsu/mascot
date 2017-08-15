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
int dimX = 2;
int dimU = 2;

/* data types for the ode solver */
typedef std::array<double,2> X_type;
typedef std::array<double,2> U_type;

/* sampling time */
double tau = 1;
/* number of intermediate steps in the ode solver */
int nint = 5;
OdeSolver ode_solver(dimX,nint,tau);

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u) -> void {

//    /* the ode describing the simple */
//    auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
//        xx[0] = u[0];
//    };
//    ode_solver(rhs,x,u);
    x[0] += u[0];
    x[1] += u[1];

};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u) -> void {
    r[0] = 0.6;
    r[1] = 0.6;
};

// forward declarations
SymbolicSet simpleCreateStateSpace(Cudd &mgr);
SymbolicSet simpleCreateInputSpace(Cudd &mgr);

auto simpleAddObstacles = [](SymbolicSet* X) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double h1[4] = {5, 4, 2, 5};
    X->remPolytope(4, H, h1, OUTER);
    double h2[4] = {10, -6, 2, 5};
    X->remPolytope(4, H, h2, OUTER);
    double h3[4] = {2, 8, 7, -4};
    X->remPolytope(4, H, h3, OUTER);
};

auto simpleConstructGoal = [](SymbolicSet* G) -> void {
    double H[2*2] = {0.5, 0,
                     0, 0.5};
    double c[2] = {-7, -7};
    G->addEllipsoid(H, c, INNER);
};

int main() {
    TicToc tt;
    double lbX[dimX]={-10, -10};
    double ubX[dimX]={10, 10};
    double etaX[dimX]= {0.6, 0.6};
    double lbU[dimU]={-1, -1};
    double ubU[dimU]={1, 1};
    double etaU[dimU]={1, 1};

    double q[dimX] = {7, 8};

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX,
                                 dimU, lbU, ubU, etaU,
                                 q);
    tt.tic();
    abs.reach(sysNext, radNext, simpleAddObstacles, simpleConstructGoal);
    tt.toc();
}
