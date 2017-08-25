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
#define a0DimX 3
#define a1DimX 3
#define a2DimX 3
#define a3DimX 3
#define aDimU 1


const int numAbs = 2;
const double tau = 0.9;
const double tauRatio = 3;
const int nSubInt = 5;

const double bw[bDimX] = {0, 0};

/* data types for the ode solver */
typedef std::array<double, bDimX> bX_t;
typedef std::array<double, bDimU> bU_t;
typedef std::array<double, a0DimX> a0X_t;
typedef std::array<double, a1DimX> a1X_t;
typedef std::array<double, a2DimX> a2X_t;
typedef std::array<double, a3DimX> a3X_t;
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
        da0xdt[0] = 0;
        da0xdt[1] = a0x[2];
        da0xdt[2] = -(a0x[1]-3.5);
    };
    solver(ODE, a0x, au);
};

auto aux0RadNext= [](a0X_t &a0r, aU_t &au, double tau, OdeSolver solver) -> void {
    a0r[0] = 0;
    a0r[1] = 0;
    a0r[2] = 0;
};

auto aux1SysNext = [](a1X_t &a1x, aU_t &au, double tau, OdeSolver solver) -> void {
    auto ODE = [](a1X_t &da1xdt, const a1X_t &a1x, aU_t &au) -> void {
        da1xdt[0] = 0;
        da1xdt[1] = a1x[2];
        da1xdt[2] = -(a1x[1]-(-3.5));
    };
    solver(ODE, a1x, au);
};

auto aux1RadNext= [](a1X_t &a1r, aU_t &au, double tau, OdeSolver solver) -> void {
    a1r[0] = 0;
    a1r[1] = 0;
    a1r[2] = 0;
};

auto aux2SysNext = [](a2X_t &a2x, aU_t &au, double tau, OdeSolver solver) -> void {
    auto ODE = [](a2X_t &da2xdt, const a2X_t &a2x, aU_t &au) -> void {
        da2xdt[0] = a2x[2];
        da2xdt[1] = 0;
        da2xdt[2] = -(a2x[0]-3.5);
    };
    solver(ODE, a2x, au);
};

auto aux2RadNext= [](a2X_t &a2r, aU_t &au, double tau, OdeSolver solver) -> void {
    a2r[0] = 0;
    a2r[1] = 0;
    a2r[2] = 0;
};

auto aux3SysNext = [](a3X_t &a3x, aU_t &au, double tau, OdeSolver solver) -> void {
    auto ODE = [](a3X_t &da3xdt, const a3X_t &a3x, aU_t &au) -> void {
        da3xdt[0] = a3x[2];
        da3xdt[1] = 0;
        da3xdt[2] = -(a3x[0]+3.5);
    };
    solver(ODE, a3x, au);
};

auto aux3RadNext= [](a3X_t &a3r, aU_t &au, double tau, OdeSolver solver) -> void {
    a3r[0] = 0;
    a3r[1] = 0;
    a3r[2] = 0;
};

auto baseAddO = [](SymbolicSet* O) -> void {
    ;
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

auto prod2AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba2x)->bool {
        return sqrt(pow(ba2x[0]-ba2x[2], 2) + pow(ba2x[1]-ba2x[3], 2)) <= 0.5;
    };
    G->addByFunction(f);
};

auto prod3AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba3x)->bool {
        return sqrt(pow(ba3x[0]-ba3x[2], 2) + pow(ba3x[1]-ba3x[3], 2)) <= 0.5;
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

    double bEtaRatioX[bDimX] = {3, 3};

    aU_t au;

    double a0LbX[a0DimX] = {   0,   2, -1.5};
    double a0UbX[a0DimX] = { 0.2,   5,  1.5};
    double a0EtaX[a0DimX] = {0.2, 0.4,  0.4};
    a0X_t a0x;
    double a0EtaRatioX[a0DimX] = {1, 2, 2};

    double a1LbX[a1DimX] = {-0.2,  -5, -1.5};
    double a1UbX[a1DimX] = {   0,  -2,  1.5};
    double a1EtaX[a1DimX] = {0.2, 0.4,  0.4};
    a1X_t a1x;
    double a1EtaRatioX[a1DimX] = {1, 2, 2};

    double a2LbX[a2DimX] = {   2, -0.2, -1.5};
    double a2UbX[a2DimX] = {   5,    0,  1.5};
    double a2EtaX[a2DimX] = {0.4,  0.2,  0.4};
    a2X_t a2x;
    double a2EtaRatioX[a2DimX] = {2, 1, 2};

    double a3LbX[a3DimX] = {  -5,    0, -1.5};
    double a3UbX[a3DimX] = {  -2,  0.2,  1.5};
    double a3EtaX[a3DimX] = {0.4,  0.2,  0.4};
    a3X_t a3x;
    double a3EtaRatioX[a3DimX] = {2, 1, 2};

    System base(bDimX, bLbX, bUbX, bEtaX, tau,
               bDimU, bLbU, bUbU, bEtaU,
               bEtaRatioX, tauRatio, nSubInt, numAbs);

    System aux0(a0DimX, a0LbX, a0UbX, a0EtaX, tau,
                a0EtaRatioX);
    System aux1(a1DimX, a1LbX, a1UbX, a1EtaX, tau,
                a1EtaRatioX);
    System aux2(a2DimX, a2LbX, a2UbX, a2EtaX, tau,
                a2EtaRatioX);
    System aux3(a3DimX, a3LbX, a3UbX, a3EtaX, tau,
                a3EtaRatioX);
    vector<System*> auxs;
    auxs.push_back(&aux0);
    auxs.push_back(&aux1);
    auxs.push_back(&aux2);
    auxs.push_back(&aux3);

    GBuchi abs("whirlpool4.txt");

    abs.initialize(&base, auxs);

    abs.initializeProdGs(prod0AddG, 0);
    abs.initializeProdGs(prod1AddG, 1);
    abs.initializeProdGs(prod2AddG, 2);
    abs.initializeProdGs(prod3AddG, 3);


    abs.initializeBaseOs(baseAddO);

    abs.computeBaseAbstractions(baseSysNext, baseRadNext, bx, bu);
    abs.computeAuxAbstractions(aux0SysNext, aux0RadNext, a0x, au, 0);
    abs.computeAuxAbstractions(aux1SysNext, aux1RadNext, a1x, au, 1);
    abs.computeAuxAbstractions(aux2SysNext, aux2RadNext, a2x, au, 2);
    abs.computeAuxAbstractions(aux3SysNext, aux3RadNext, a3x, au, 3);
    abs.composeAbstractions();

    int startAbs = 1;
    int minToGoCoarser = 2;
    int minToBeValid = 5;

    abs.gBuchi(startAbs, minToGoCoarser, minToBeValid);
}


int main() {
    composition();
}
