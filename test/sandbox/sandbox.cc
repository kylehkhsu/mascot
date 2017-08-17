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
#define dimB 3

const double w[dimX + dimB] = {0, 0, 0, 0, 0};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;
typedef std::array<double, dimB> B_type;
typedef std::array<double, dimX + dimB> XB_type;
typedef std::array<double, 1> Udum_type;

auto sysNextTog = [](XB_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](XB_type &xx, const XB_type &x, U_type &u) -> void {
        xx[0] = (-0.5 * x[0] + 1 * x[1]) * u[1];
        xx[1] = (-0.5 * x[0] + 0.5 * x[1] + u[0]) * u[1];
        xx[2] = 0;
        xx[3] = x[4];
        xx[4] = -(x[3]-3.5);
    };
    solver(ODE, x, u);
};

auto radNextTog = [](XB_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](XB_type &drdt, const XB_type &r, const U_type &u) -> void {
        drdt[0] = -0.25*r[0] + 1*r[1] + w[0];
        drdt[1] = -0.25*r[0] + 0.5*r[1] + w[1];
        drdt[2] = 0;
        drdt[3] = 0;
        drdt[4] = 0;
    };
    solver(radODE, r, u);
    r[2] = 0;
    r[3] = 0;
    r[4] = 0;
};

auto sysNextSepX = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](X_type &xx, const X_type &x, U_type &u) -> void {
        xx[0] = (-0.5 * x[0] + 1 * x[1]) * u[1];
        xx[1] = (-0.5 * x[0] + 0.5 * x[1] + u[0]) * u[1];
    };
    solver(ODE, x, u);
};

auto radNextSepX = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        drdt[0] = -0.25*r[0] + 1*r[1] + w[0];
        drdt[1] = -0.25*r[0] + 0.5*r[1] + w[1];
    };
    solver(radODE, r, u);
};

auto sysNextSepB = [](B_type &x, Udum_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](B_type &xx, const B_type &x, Udum_type &u) -> void {
        xx[0] = 0;
        xx[1] = x[2];
        xx[2] = -(x[1]-3.5);
    };
    solver(ODE, x, u);
};

auto radNextSepB = [](B_type &r, Udum_type &u, double tau, OdeSolver solver) -> void {
    r[0] = 0;
    r[1] = 0;
    r[2] = 0;
};

void together() {
    TicToc tt;
    tt.tic();
    double lbX[dimX + dimB]={-6, -6, -0.1, 2, -1.5};
    double ubX[dimX + dimB]={ 6,  6,  0.1, 5,  1.5};
    double etaX[dimX + dimB]= {0.6, 0.6, 0.2, 0.2, 0.2};

    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2,   1};
    double etaU[dimU]= {0.5, 0.2};

    double tau = 0.9;

    Cudd mgr;

    OdeSolver solver(dimX + dimB, 5, tau);

    SymbolicSet X(mgr, dimX + dimB, lbX, ubX, etaX, tau);
    X.addGridPoints();
//    X.printInfo(1);
    SymbolicSet U(mgr, dimU, lbU, ubU, etaU, 0);
    U.addGridPoints();
//    U.printInfo(1);
    SymbolicSet X2(X, 1);

    SymbolicModelGrowthBound<XB_type, U_type> Abs(&X, &U, &X2);

    Abs.computeTransitionRelation(sysNextTog, radNextTog, solver);

    SymbolicSet T = Abs.getTransitionRelation();
    tt.toc();
    T.printInfo(1);
}

void separate() {
    TicToc tt;
    tt.tic();
    double lbX[dimX]={-6, -6};
    double ubX[dimX]={ 6,  6 };
    double etaX[dimX]= {0.6, 0.6};
    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2,   1};
    double etaU[dimU]= {0.5, 0.2};

    double tau = 0.9;

    Cudd mgr;
    OdeSolver solver(dimX, 5, tau);

    SymbolicSet X(mgr, dimX, lbX, ubX, etaX, tau);
    X.addGridPoints();
    SymbolicSet U(mgr, dimU, lbU, ubU, etaU, tau);
    U.addGridPoints();
    SymbolicSet X2(X, 1);

    SymbolicModelGrowthBound<X_type, U_type> Abs(&X, &U, &X2);
    Abs.computeTransitionRelation(sysNextSepX, radNextSepX, solver);
    SymbolicSet T = Abs.getTransitionRelation();
    T.printInfo(1);

    double lbB[dimB] = {-0.1, 2, -1.5};
    double ubB[dimB] = { 0.1, 5,  1.5};
    double etaB[dimB] = {0.2, 0.2, 0.2};
    double lbUdum[1] = {-0.5};
    double ubUdum[1] = {0.5};
    double etaUdum[1] = {1};

    SymbolicSet B(mgr, dimB, lbB, ubB, etaB, tau);
    B.addGridPoints();
    SymbolicSet Udum(mgr, 1, lbUdum, ubUdum, etaUdum, 0);
    Udum.addGridPoints();
//    Udum.printInfo(1);
    SymbolicSet B2(B, 1);

    SymbolicModelGrowthBound<B_type, Udum_type> ball(&B, &Udum, &B2);
    ball.computeTransitionRelation(sysNextSepB, radNextSepB, solver);
    SymbolicSet TB = ball.getTransitionRelation();
    TB.printInfo(1);

    SymbolicSet All(T, TB);
    All.addGridPoints();
    All.symbolicSet_ &= T.symbolicSet_ & TB.symbolicSet_;
    tt.toc();
    All.printInfo(1);

}



int main() {
//    together();
    separate();






}
