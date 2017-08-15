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
    x[0] += u[0];
    x[1] += u[1];
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u) -> void {
    r[0] = r[0] + abs(u[0]) * 0.4;
    r[1] = r[1] + abs(u[1]) * 0.4;
};

// forward declarations
SymbolicSet simpleCreateStateSpace(Cudd &mgr);
SymbolicSet simpleCreateInputSpace(Cudd &mgr);

auto simpleAddObstacles = [](SymbolicSet* X) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double h1[4] = {-25, 27, -5, 30};
    X->remPolytope(4, H, h1, OUTER);
    double h2[4] = {-5, 28.5, -5, 8};
    X->remPolytope(4, H, h2, OUTER);
    double h3[4] = {-5, 8, -8, 25};
    X->remPolytope(4, H, h3, OUTER);
    double h4[4] = {-8, 20, -23, 25};
    X->remPolytope(4, H, h4, OUTER);
};

auto simpleConstructGoal = [](SymbolicSet* G) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {-10, 15, -15, 20};
    G->addPolytope(4, H, c, INNER);
};

int main() {
    TicToc tt;
    Cudd mgr;
    double lbX[dimX]={0, 0};
    double ubX[dimX]={30, 30};
    double etaX[dimX]= {1, 1};
    double lbU[dimU]={-2, -2};
    double ubU[dimU]={2, 2};
    double etaU[dimU]= {0.5, 0.5};

    double q[dimX] = {29, 29};

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX,
                                 dimU, lbU, ubU, etaU,
                                 q);
    tt.tic();
    abs.reach(sysNext, radNext, simpleAddObstacles, simpleConstructGoal);
    tt.toc();
}
