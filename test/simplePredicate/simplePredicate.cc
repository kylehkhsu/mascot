#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "Reach.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"
#include "Product.hh"

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


auto dynSysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](X_type &xx, const X_type &x, U_type &u) -> void {
        xx[0] = (-0.5 * x[0] + 1 * x[1]) * u[1];
        xx[1] = (-0.5 * x[0] + 0.5 * x[1] + u[0]) * u[1];
    };
    solver(ODE, x, u);
};

auto dynRadNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        drdt[0] = -0.25*r[0] + 1*r[1] + w[0];
        drdt[1] = -0.25*r[0] + 0.5*r[1] + w[1];
    };
    solver(radODE, r, u);
};

auto ballSysNext = [](B_type &x, Y_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](B_type &xx, const B_type &x, Y_type &u) -> void {
        xx[0] = 0;
        xx[1] = x[2];
        xx[2] = -(x[1]-3.5);
    };
    solver(ODE, x, u);
};

auto ballRadNext= [](B_type &r, Y_type &u, double tau, OdeSolver solver) -> void {
    r[0] = 0;
    r[1] = 0;
    r[2] = 0;
};

auto dummyAddO = [](SymbolicSet* O) -> void {
    ;
};

auto whirlpoolAddG = [](SymbolicSet* G) -> void {
    auto f = [](double* x)->bool {
        return sqrt(pow(x[0]-x[2], 2) + pow(x[1]-x[3], 2)) <= 1;
    };
    G->addByFunction(f);
};

auto whirlpoolAddI = [](SymbolicSet* I) -> void {
    double q[5] = {1, 1, 0.1, 4.6, 0};
    I->addPoint(q);
};


void product() {
    double lbX[dimX]={-6, -6};
    double ubX[dimX]={ 6,  6 };
    double etaX[dimX]= {0.6, 0.6};
    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2,   1};
    double etaU[dimU]= {0.5, 0.2};
    X_type x;
    U_type u;

    double dynEtaRatio[dimX] = {3, 3};
    double dynTauRatio = 1;

    double lbB[dimB] = {   0,   2, -1.5};
    double ubB[dimB] = { 0.2,   5,  1.5};
    double etaB[dimB] = {0.2, 0.2,  0.2};
    double lbY[1] = { 0};
    double ubY[1] = { 1};
    double etaY[1] = {1};
    B_type b;
    Y_type y;

    double ballEtaRatio[dimB] = {1, 1, 1};
    double ballTauRatio = 1;

    System dyn(dimX, lbX, ubX, etaX, tau,
               dimU, lbU, ubU, etaU,
               dynEtaRatio, dynTauRatio, nint, numAbs);

    System ball(dimB, lbB, ubB, etaB, tau,
                dimY, lbY, ubY, etaY,
                ballEtaRatio, ballTauRatio, nint, numAbs);
    vector<System*> balls;
    balls.push_back(&ball);

    Reach<X_type, U_type> abs("product.txt");

    abs.initializeProduct(&dyn, balls);
    abs.computeDynAbstractions(dynSysNext, dynRadNext, x, u);
    abs.computePredAbstractions(ballSysNext, ballRadNext, b, y, 0);
    abs.computeProducts();

    int readXX = 0;
    int readAbs = 1; // don't change
    abs.initializeAdaptive(readXX, readAbs, dummyAddO);
    abs.initializeReach(whirlpoolAddG, whirlpoolAddI);
    abs.computeAbstractions(dynSysNext, dynRadNext);

    int startAbs = 0;
    int minToGoCoarser = 2;
    int minToBeValid = 5;
    int earlyBreak = 1;

    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);

}


int main() {
    product();
//    separate();




}
