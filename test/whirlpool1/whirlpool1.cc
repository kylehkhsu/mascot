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
#define aDimU 1

const int numAbs = 2;
const double tau = 0.5;
const double tauRatio = 2;
const int nSubInt = 5;

const double bw[bDimX] = {0, 0};

/* data types for the ode solver */
typedef std::array<double, bDimX> bX_t;
typedef std::array<double, bDimU> bU_t;
typedef std::array<double, a0DimX> a0X_t;
typedef std::array<double, 1> aU_t;

auto baseSysNext = [](bX_t &bx, bU_t &bu, double tau, OdeSolver solver) -> void {
    auto ODE = [](bX_t &dbxdt, const bX_t &bx, bU_t &bu) -> void {
        dbxdt[0] = (-0.5 * bx[0] + 1 * bx[1]) * bu[1];
        dbxdt[1] = (-0.5 * bx[0] + 0.5 * bx[1] + bu[0]) * bu[1];
    };
    solver(ODE, bx, bu);
};

auto baseRadNext = [](bX_t &br, bU_t &bu, double tau, OdeSolver solver) -> void {
    auto radODE = [](bX_t &dbrdt, const bX_t &br, const bU_t &bu) -> void {
        dbrdt[0] = -0.25*br[0] + 1*br[1] + bw[0];
        dbrdt[1] = -0.25*br[0] + 0.5*br[1] + bw[1];
    };
    solver(radODE, br, bu);
};

auto aux0SysNext = [](a0X_t &a0x, aU_t &au, double tau, OdeSolver solver) -> void {
    auto ODE = [](a0X_t &da0xdt, const a0X_t &a0x, aU_t &au) -> void {
        da0xdt[0] = a0x[1];
        da0xdt[1] = -(a0x[0]-3.5);
    };
    solver(ODE, a0x, au);
};

auto aux0RadNext= [](a0X_t &a0r, aU_t &au, double tau, OdeSolver solver) -> void {
//    auto ODE = [](a0X_t &da0rdt, const a0X_t &a0r, aU_t &au) -> void {
//        da0rdt[0] = a0r[1];
//        da0rdt[1] = a0r[0];
//    };
//    solver(ODE, a0r, au);
    a0r[0] = 0;
    a0r[1] = 0;
};
auto baseAddO = [](SymbolicSet* O) -> void {
    ;
};

auto prod0AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba0x)->bool {
        return sqrt(pow(ba0x[0]-0, 2) + pow(ba0x[1]-ba0x[2], 2)) <= 0.5;
    };
    G->addByFunction(f);
};

void composition() {
    double bLbX[bDimX]={-6, -6};
    double bUbX[bDimX]={ 6,  6 };
    double bEtaX[bDimX]= {0.6, 0.6};
    double bLbU[bDimU]= {-2, 0.5};
    double bUbU[bDimU]= { 2,   1};
    double bEtaU[bDimU]= {0.5, 0.2};
    bX_t bx;
    bU_t bu;

    double bEtaRatioX[bDimX] = {2, 2};

    double a0LbX[a0DimX] = {   0, -1.5};
    double a0UbX[a0DimX] = {   7,  1.5};
    double a0EtaX[a0DimX] = {0.4,  0.4};

    a0X_t a0x;
    aU_t au;

    double a0EtaRatioX[a0DimX] = {2, 2};

    System base(bDimX, bLbX, bUbX, bEtaX, tau,
               bDimU, bLbU, bUbU, bEtaU,
               bEtaRatioX, tauRatio, nSubInt, numAbs);

    System aux0(a0DimX, a0LbX, a0UbX, a0EtaX, tau,
                a0EtaRatioX);
    vector<System*> auxs;
    auxs.push_back(&aux0);

    GBuchi abs("whirlpool1.txt");

    abs.initialize(&base, auxs);

    abs.initializeProdGs(prod0AddG, 0);
    abs.initializeBaseOs(baseAddO);

    abs.computeBaseAbstractions(baseSysNext, baseRadNext, bx, bu);
    abs.computeAuxAbstractions(aux0SysNext, aux0RadNext, a0x, au, 0);
    abs.composeAbstractions();

    int startAbs = 1;
    int minToGoCoarser = 2;
    int minToBeValid = 5;

    abs.gBuchi(startAbs, minToGoCoarser, minToBeValid);
}


int main() {
    composition();
}
