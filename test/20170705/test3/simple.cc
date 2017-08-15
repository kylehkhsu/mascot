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
double tau = 0.5;
/* number of intermediate steps in the ode solver */
int nint = 5;
OdeSolver ode_solver(dimX,nint,tau);

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [tau](X_type &x, U_type &u) -> void {
    x[0] += u[0] * tau;
    x[1] += u[1] * tau;
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [tau](X_type &r, U_type &u) -> void {
    r[0] = r[0] + abs(u[0]) * tau / 2;
    r[1] = r[1] + abs(u[1]) * tau / 2;
};

auto radNo = [](X_type &r, U_type &u) -> void {
    r[0] = 0;
    r[1] = 0;
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
    double h2[4] = {-5, 27, -5, 8};
    X->remPolytope(4, H, h2, OUTER);
    double h3[4] = {-5, 8, -8, 25};
    X->remPolytope(4, H, h3, OUTER);
    double h4[4] = {-8, 20, -23, 25};
    X->remPolytope(4, H, h4, OUTER);
    double h5[4] = {-16, 18, -12, 24};
    X->remPolytope(4, H, h5, OUTER);
    double h6[4] = {-11, 17, -12, 14};
    X->remPolytope(4, H, h6, OUTER);
};

auto simpleConstructGoal = [](SymbolicSet* G) -> void {
    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {-10, 15, -15, 20};
    G->addPolytope(4, H, c, INNER);
};

auto simpleConstructInit = [](SymbolicSet* I) -> void {
    double q[2] = {23, 27};
    I->addPoint(q);
};

int main() {
    TicToc tt;
    Cudd mgr;
    double lbX[dimX]={0, 0};
    double ubX[dimX]={40, 40};
    double etaX[dimX]= {0.44444, 0.44444};
    double lbU[dimU]={-2, -2};
    double ubU[dimU]={2, 2};
    double etaU[dimU]= {0.5, 0.5};

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX,
                                 dimU, lbU, ubU, etaU);
    tt.tic();
    abs.reach(sysNext, radNext, simpleAddObstacles, simpleConstructGoal, simpleConstructInit);
    tt.toc();
//    tt.tic();
//        abs.preProcess(sysNext, radNext, radNo, simpleAddObstacles, simpleConstructGoal, simpleConstructInit);
//    tt.toc();
}
