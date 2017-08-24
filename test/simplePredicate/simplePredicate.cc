#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"
#include "GBuchi.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define dimX 2
#define dimU 2
#define dimB 3
#define dimY 1

const int numAbs = 2;
const double tau = 0.9;
const double tauRatio = 2;
const int nSubInt = 5;

const double w[dimX + dimB] = {0, 0, 0, 0, 0};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;
typedef std::array<double, dimB> B_type;
typedef std::array<double, dimX + dimB> XB_type;
typedef std::array<double, 1> Y_type;


auto baseSysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](X_type &xx, const X_type &x, U_type &u) -> void {
        xx[0] = (-0.5 * x[0] + 1 * x[1]) * u[1];
        xx[1] = (-0.5 * x[0] + 0.5 * x[1] + u[0]) * u[1];
    };
    solver(ODE, x, u);
};

auto baseRadNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        drdt[0] = -0.25*r[0] + 1*r[1] + w[0];
        drdt[1] = -0.25*r[0] + 0.5*r[1] + w[1];
    };
    solver(radODE, r, u);
};

auto aux0SysNext = [](B_type &x, Y_type &u, double tau, OdeSolver solver) -> void {
    auto ODE = [](B_type &xx, const B_type &x, Y_type &u) -> void {
        xx[0] = 0;
        xx[1] = x[2];
        xx[2] = -(x[1]-3.5);
    };
    solver(ODE, x, u);
};

auto aux0RadNext= [](B_type &r, Y_type &u, double tau, OdeSolver solver) -> void {
    r[0] = 0;
    r[1] = 0;
    r[2] = 0;
};

auto baseAddO = [](SymbolicSet* O) -> void {
    ;
};

auto prod0AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* x)->bool {
        return sqrt(pow(x[0]-x[2], 2) + pow(x[1]-x[3], 2)) <= 0.5;
    };
    G->addByFunction(f);
};

auto baseAddI = [](SymbolicSet* I) -> void {
    double q[dimX] = {1, 1};
    I->addPoint(q);
};

auto aux0AddI = [](SymbolicSet* I) -> void {
    double q[dimB] = {0.1, 4, 0};
    I->addPoint(q);
};


void composition() {
    double lbX[dimX]={-6, -6};
    double ubX[dimX]={ 6,  6 };
    double etaX[dimX]= {0.6, 0.6};
    double lbU[dimU]= {-2, 0.5};
    double ubU[dimU]= { 2,   1};
    double etaU[dimU]= {0.5, 0.2};
    X_type x;
    U_type u;

    double baseEtaRatio[dimX] = {2, 2};

    double lbB[dimB] = {   0,   2, -1.5};
    double ubB[dimB] = { 0.2,   5,  1.5};
    double etaB[dimB] = {0.2, 0.2,  0.2};

    B_type b;
    Y_type y;

    double ballEtaRatio[dimB] = {1, 1, 1};


    System base(dimX, lbX, ubX, etaX, tau,
               dimU, lbU, ubU, etaU,
               baseEtaRatio, tauRatio, nSubInt, numAbs);

    System ball(dimB, lbB, ubB, etaB, tau,
                ballEtaRatio);
    vector<System*> balls;
    balls.push_back(&ball);

    GBuchi abs("composition.txt");

    abs.initialize(&base, balls);

    abs.initializeProdGs(prod0AddG, 0);
    abs.initializeBaseIs(baseAddI);
    abs.initializeAuxIs(aux0AddI, 0);
    abs.initializeBaseOs(baseAddO);

    abs.computeBaseAbstractions(baseSysNext, baseRadNext, x, u);
    abs.computeAuxAbstractions(aux0SysNext, aux0RadNext, b, y, 0);
    abs.composeAbstractions();

    int startAbs = 1;
    int minToGoCoarser = 2;
    int minToBeValid = 5;
//    int earlyBreak = 1;

//    abs.reachOne(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
    abs.alwaysEventuallyOne(startAbs, minToGoCoarser, minToBeValid);
}


int main() {
    composition();
//    separate();




}
