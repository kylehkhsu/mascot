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

#define bDimX 2
#define bDimU 2
#define a0DimX 2
#define a1DimX 2
#define aDimU 1


const int numAbs = 2;
const double tau = 0.6;
const double tauRatio = 2;
const int nSubInt = 5;

const double bw[bDimX] = {0.05, 0.05};

/* data types for the ode solver */
typedef std::array<double, bDimX> bX_t;
typedef std::array<double, bDimU> bU_t;
typedef std::array<double, a0DimX> a0X_t;
typedef std::array<double, a1DimX> a1X_t;
typedef std::array<double, 1> aU_t;

auto baseSysNext = [](bX_t &bx, bU_t &bu, double tau, OdeSolver solver) -> void {
    auto ODE = [](bX_t &dbxdt, const bX_t &bx, bU_t &bu) -> void {
        dbxdt[0] = bu[0];
        dbxdt[1] = bu[1];
    };
    solver(ODE, bx, bu);
};

auto baseRadNext = [](bX_t &br, bU_t &bu, double tau, OdeSolver solver) -> void {
    auto radODE = [](bX_t &dbrdt, const bX_t &br, const bU_t &bu) -> void {
        dbrdt[0] = bw[0];
        dbrdt[1] = bw[1];
    };
    solver(radODE, br, bu);
};

auto aux0SysNext = [](a0X_t &a0x, aU_t &au, double tau, OdeSolver solver) -> void {
    a0x[0] = a0x[0] + 0.5;
    if (a0x[0] == 3.25) {
        a0x[0] = -2.75;
    }
};

auto aux0RadNext= [](a0X_t &a0r, aU_t &au, double tau, OdeSolver solver) -> void {
    a0r[0] = 0;
    a0r[1] = 0;
};

auto aux1SysNext = [](a1X_t &a1x, aU_t &au, double tau, OdeSolver solver) -> void {
    a1x[0] = a1x[0] + 1;
    if (a1x[0] == 3.5) {
        a1x[0] = -2.5;
    }
};

auto aux1RadNext= [](a1X_t &a1r, aU_t &au, double tau, OdeSolver solver) -> void {
    a1r[0] = 0;
    a1r[1] = 0;
};

auto baseAddO = [](SymbolicSet* O) -> void {
    double H[4*2]={-1, 0,
                    1, 0,
                    0,-1,
                    0, 1};
    double h1[4] = {7, -2.8, 1, 1};
    double h2[4] = {1.5, 1.5, 1, 1};
    double h3[4] = {-2.8, 7, 1, 1};
    O->addPolytope(4, H, h1, OUTER);
    O->addPolytope(4, H, h2, OUTER);
    O->addPolytope(4, H, h3, OUTER);
};

auto prod0AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba0x)->bool {
        return sqrt(pow(ba0x[0]-ba0x[2], 2) + pow(ba0x[1]-ba0x[3], 2)) <= 0.5;
    };
    G->addByFunction(f);
};

auto prod1AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba1x)->bool {
        return sqrt(pow(ba1x[0]-ba1x[2], 2) + pow(ba1x[1]-ba1x[3], 2)) <= 0.5;
    };
    G->addByFunction(f);
};

void composition() {
    double bLbX[bDimX]={-6,  -6 };
    double bUbX[bDimX]={ 6,   6 };
    double bEtaX[bDimX]= {0.6, 0.6};
    double bLbU[bDimU]= {-1.3, -1.3};
    double bUbU[bDimU]= { 1.3,  1.3};
    double bEtaU[bDimU]= {0.5, 0.5};
    bX_t bx;
    bU_t bu;

    double bEtaRatioX[bDimX] = {2, 2};

    aU_t au;

    double a0LbX[a0DimX] = {  -3, 0};
    double a0UbX[a0DimX] = {   3, 8};
    double a0EtaX[a0DimX] = { 0.5, 8};
    a0X_t a0x;
    double a0EtaRatioX[a0DimX] = {1, 1};

    double a1LbX[a1DimX] = {  -3, -8};
    double a1UbX[a1DimX] = {   3,  0};
    double a1EtaX[a1DimX] = {  1,  8};
    a1X_t a1x;
    double a1EtaRatioX[a1DimX] = {1, 1};

    System base(bDimX, bLbX, bUbX, bEtaX, tau,
               bDimU, bLbU, bUbU, bEtaU,
               bEtaRatioX, tauRatio, nSubInt, numAbs);

    System aux0(a0DimX, a0LbX, a0UbX, a0EtaX, tau,
                a0EtaRatioX);
    System aux1(a1DimX, a1LbX, a1UbX, a1EtaX, tau,
                a1EtaRatioX);
    vector<System*> auxs;
    auxs.push_back(&aux0);
    auxs.push_back(&aux1);

    GBuchi abs("simpleBall2.txt");

    abs.initialize(&base, auxs);

    abs.initializeProdGs(prod0AddG, 0);
    abs.initializeProdGs(prod1AddG, 1);

    abs.initializeBaseOs(baseAddO);

    abs.computeBaseAbstractions(baseSysNext, baseRadNext, bx, bu);
    abs.computeAuxAbstractions(aux0SysNext, aux0RadNext, a0x, au, 0);
    abs.computeAuxAbstractions(aux1SysNext, aux1RadNext, a1x, au, 1);
    abs.composeAbstractions();

    int startAbs = 1;
    int minToGoCoarser = 3;
    int minToBeValid = 3;

    abs.gBuchi(startAbs, minToGoCoarser, minToBeValid);
}


int main() {
    composition();
}
