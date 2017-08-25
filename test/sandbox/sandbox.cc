#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include <tuple>

#include "cuddObj.hh"

//#include "Adaptive.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"


using namespace std;
using namespace scots;
using namespace helper;

#define dimX 2
#define dimU 2
#define dimB 3
#define dimY 1

const int numAbs = 1;
const double tau = 0.9;
const int nint = 5;

const double w[dimX + dimB] = {0, 0, 0, 0, 0};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;
typedef std::array<double, dimB> B_type;
typedef std::array<double, dimX + dimB> XB_type;
typedef std::array<double, 1> Y_type;

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

auto sysNextSepB = [](B_type &x, Y_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](B_type &xx, const B_type &x, Y_type &u) -> void {
        xx[0] = 0;
        xx[1] = x[2];
        xx[2] = -(x[1]-3.5);
    };
    solver(ODE, x, u);
};

auto radNextSepB = [](B_type &r, Y_type &u, double tau, OdeSolver solver) -> void {
    r[0] = 0;
    r[1] = 0;
    r[2] = 0;
};

auto dummyAddO = [](SymbolicSet* O) -> void {
    ;
};

void together() {
//    TicToc tt;
//    tt.tic();
//    double lbX[dimX + dimB]={-6, -6,    0, 2, -1.5};
//    double ubX[dimX + dimB]={ 6,  6,  0.2, 5,  1.5};
//    double etaX[dimX + dimB]= {0.6, 0.6, 0.2, 0.2, 0.2};

//    double lbU[dimU]= {-2, 0.5};
//    double ubU[dimU]= { 2,   1};
//    double etaU[dimU]= {0.5, 0.2};

//    double tau = 0.9;

//    Cudd mgr;

//    OdeSolver solver(dimX + dimB, 5, tau);

//    SymbolicSet X(mgr, dimX + dimB, lbX, ubX, etaX, tau);
//    X.addGridPoints();
//    X.printInfo(1);
//    SymbolicSet X2(X, 1);
//    SymbolicSet U(mgr, dimU, lbU, ubU, etaU, 0);
//    U.addGridPoints();
//    X2.printInfo(1);
//    U.printInfo(1);


//    SymbolicModelGrowthBound<XB_type, U_type> Abs(&X, &U, &X2);

//    Abs.computeTransitionRelation(sysNextTog, radNextTog, solver);

//    SymbolicSet T = Abs.getTransitionRelation();
//    tt.toc();
//    T.printInfo(1);
}

void separate() {
//    Cudd mgr;
//    double lbX[dimX]={-6, -6};
//    double ubX[dimX]={ 6,  6 };
//    double etaX[dimX]= {0.6, 0.6};
//    double tau = 0.9;
//    double lbU[dimU]= {-2, 0.5};
//    double ubU[dimU]= { 2,   1};
//    double etaU[dimU]= {0.5, 0.2};

//    double lbB[dimB] = {   0,   2, -1.5};
//    double ubB[dimB] = { 0.2,   5,  1.5};
//    double etaB[dimB] = {0.2, 0.2,  0.2};
//    double lbY[1] = { 0};
//    double ubY[1] = { 1};
//    double etaY[1] = {1};

//    SymbolicSet X(mgr, dimX, lbX, ubX, etaX, tau);
//    SymbolicSet B(mgr, dimB, lbB, ubB, etaB, tau);
//    SymbolicSet X2(X, 1);
//    SymbolicSet B2(B, 1);
//    SymbolicSet U(mgr, dimU, lbU, ubU, etaU, 0);
//    SymbolicSet Y(mgr,    1, lbY, ubY, etaY, 0);

//    X.addGridPoints();
//    B.addGridPoints();
//    U.addGridPoints();
//    Y.addGridPoints();


//    OdeSolver solverD(dimX, 5, tau);
//    OdeSolver solverP(dimB, 5, tau);

//    SymbolicModelGrowthBound<X_type, U_type> absD(&X, &U, &X2);
//    SymbolicModelGrowthBound<B_type, Y_type> absP(&B, &Y, &B2);

//    absD.computeTransitionRelation(sysNextSepX, radNextSepX, solverD);
//    absP.computeTransitionRelation(sysNextSepB, radNextSepB, solverP);

//    SymbolicSet TD = absD.getTransitionRelation();
//    TD.printInfo(1);
//    SymbolicSet TP = absP.getTransitionRelation();
//    TP.printInfo(1);

//    // get 'rid' of dummy input
//    SymbolicSet TPP(B, B2);
//    TPP.symbolicSet_ = TP.symbolicSet_.ExistAbstract(Y.getCube());
//    TPP.printInfo(1);

//    // construct composed T SymbolicSet
//    SymbolicSet d1(X, B);
//    SymbolicSet d2(d1, U);
//    SymbolicSet d3(d2, X2);
//    SymbolicSet T(d3, B2);


//    T.symbolicSet_ = TD.symbolicSet_ & TPP.symbolicSet_;
//    T.printInfo(1);


}

void testRK() {
    auto sysNextSepX = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
        auto ODE = [](X_type &xx, const X_type &x, U_type &u) -> void {
            xx[0] = (-0.5 * x[0] + 1 * x[1]) * u[1];
            xx[1] = (-0.5 * x[0] + 0.5 * x[1] + u[0]) * u[1];
        };
        solver(ODE, x, u);
    };

    OdeSolver solver(2, 5, 0.45);
    X_type x;
    x[0] = 5;
    x[1] = 2;

    cout << "x: " << x[0] << ' ' << x[1] << '\n';

    U_type u;
    u[0] = 0.5;
    u[1] = 0.7;

    sysNextSepX(x,u,0.45,solver);

    cout << "x: " << x[0] << ' ' << x[1] << '\n';

}


void product() {
//    double lbX[dimX]={-6, -6};
//    double ubX[dimX]={ 6,  6 };
//    double etaX[dimX]= {0.6, 0.6};
//    double lbU[dimU]= {-2, 0.5};
//    double ubU[dimU]= { 2,   1};
//    double etaU[dimU]= {0.5, 0.2};
//    X_type x;
//    U_type u;

//    double dynEtaRatio[dimX] = {3, 3};
//    double dynTauRatio = 1;

//    double lbB[dimB] = {   0,   2, -1.5};
//    double ubB[dimB] = { 0.2,   5,  1.5};
//    double etaB[dimB] = {0.2, 0.2,  0.2};
//    double lbY[1] = { 0};
//    double ubY[1] = { 1};
//    double etaY[1] = {1};
//    B_type b;
//    Y_type y;

//    double ballEtaRatio[dimB] = {1, 1, 1};
//    double ballTauRatio = 1;

//    System dyn(dimX, lbX, ubX, etaX, tau,
//               dimU, lbU, ubU, etaU,
//               dynEtaRatio, dynTauRatio, nint, numAbs);

//    System ball(dimB, lbB, ubB, etaB, tau,
//                dimY, lbY, ubY, etaY,
//                ballEtaRatio, ballTauRatio, nint, numAbs);
//    vector<System*> balls;
//    balls.push_back(&ball);

//    Product<XB_type, U_type> abs("product.txt");

//    abs.initializeProduct(&dyn, balls);
//    abs.computeDynAbstractions(sysNextSepX, radNextSepX, x, u);
//    abs.computePredAbstractions(sysNextSepB, radNextSepB, b, y, 0);
//    abs.computeProducts();

//    int readXX = 0;
//    int readAbs = 1;
//    abs.initializeAdaptive(readXX, readAbs, dummyAddO);
//    abs.computeAbstractions(sysNextTog, radNextTog);
}


int main() {

    testRK();
//    product();
//    together();
//    separate();




}
