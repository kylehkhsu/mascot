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
#define a2DimX 2
#define aDimU 1


const int numAbs = 3; // 2-1;
const double tau = 1.2; // 0.3 for the finest layer
const double tauRatio = 2;
const int nSubInt = 5; // ??

const double bw[bDimX] = {0.05, 0.05}; // ??

/* data types for the ode solver */
typedef std::array<double, bDimX> bX_t;
typedef std::array<double, bDimU> bU_t;
typedef std::array<double, a0DimX> a0X_t;
typedef std::array<double, a1DimX> a1X_t;
typedef std::array<double, a2DimX> a2X_t;
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
        dbrdt[0] = (-0.5*br[0] + br[1]) * bu[1] + bw[0];
        dbrdt[1] = ( 0.5*br[0] + 0.5*br[1]) * bu[1] + bw[1];
    };
    solver(radODE, br, bu);
};

auto aux0SysNext = [](a0X_t &a0x, aU_t &au, double tau, OdeSolver solver) -> void {
    a0x[0] = a0x[0] + 1 * tau/0.6;
    if (a0x[0] > 2.75) {
        a0x[0] = -2.75;
    }
};

auto aux0RadNext= [](a0X_t &a0r, aU_t &au, double tau, OdeSolver solver) -> void {
    a0r[0] = 0;
    a0r[1] = 0;
};

auto aux1SysNext = [](a1X_t &a1x, aU_t &au, double tau, OdeSolver solver) -> void {
    a1x[1] = a1x[1] - 1 * tau/0.6;
    if (a1x[1] < -2.75) {
        a1x[1] = 2.75;
    }
};

auto aux1RadNext= [](a1X_t &a1r, aU_t &au, double tau, OdeSolver solver) -> void {
    a1r[0] = 0;
    a1r[1] = 0;
};

auto aux2SysNext = [](a2X_t &a2x, aU_t &au, double tau, OdeSolver solver) -> void {
    a2x[0] = a2x[0] - 1 * tau/0.6;
    if (a2x[0] < -2.75) {
        a2x[0] = 2.75;
    }
};

auto aux2RadNext= [](a2X_t &a2r, aU_t &au, double tau, OdeSolver solver) -> void {
    a2r[0] = 0;
    a2r[1] = 0;
};

auto baseAddO = [](SymbolicSet* O) -> void {
    ;
};

auto prod0AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba0x)->bool {
        return sqrt(pow(ba0x[0]-ba0x[2], 2) + pow(ba0x[1]-ba0x[3], 2)) <= 1;
    };
    G->addByFunction(f);
};

auto prod1AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba1x)->bool {
        return sqrt(pow(ba1x[0]-ba1x[2], 2) + pow(ba1x[1]-ba1x[3], 2)) <= 1;
    };
    G->addByFunction(f);
};

auto prod2AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba2x)->bool {
        return sqrt(pow(ba2x[0]-ba2x[2], 2) + pow(ba2x[1]-ba2x[3], 2)) <= 1;
    };
    G->addByFunction(f);
};


void composition() {
    double bLbX[bDimX]={-6, -6};
    double bUbX[bDimX]={ 6,  6 };
    double bEtaX[bDimX]= {1.2, 1.2}; // 0.3 for the finest layer
    double bLbU[bDimU]= {-2, 0.5};
    double bUbU[bDimU]= { 2,   1};
    double bEtaU[bDimU]= {0.5, 0.2};
    bX_t bx;
    bU_t bu;

    double bEtaRatioX[bDimX] = {2, 2};

    aU_t au;

    double a0LbX[a0DimX] = {  -3, 0};
    double a0UbX[a0DimX] = {   3, 8};
    double a0EtaX[a0DimX] = {  1, 8}; // {0.5 8} for the finest layer
    a0X_t a0x;
    double a0EtaRatioX[a0DimX] = {2, 1};

    double a1LbX[a1DimX] = {   0, -3};
    double a1UbX[a1DimX] = {   8,  3};
    double a1EtaX[a1DimX] = {  8,  1}; // {8, 0.5} for the finest layer
    a1X_t a1x;
    double a1EtaRatioX[a1DimX] = {1, 2};

    double a2LbX[a2DimX] = {  -3, -8};
    double a2UbX[a2DimX] = {   3,  0};
    double a2EtaX[a2DimX] = {  1,  8}; // {0.5 8} for the finest layer
    a2X_t a2x;
    double a2EtaRatioX[a2DimX] = {2, 1};

    System base(bDimX, bLbX, bUbX, bEtaX, tau,
               bDimU, bLbU, bUbU, bEtaU,
               bEtaRatioX, tauRatio, nSubInt, numAbs);

    System aux0(a0DimX, a0LbX, a0UbX, a0EtaX, tau,
                a0EtaRatioX);
    System aux1(a1DimX, a1LbX, a1UbX, a1EtaX, tau,
                a1EtaRatioX);
    System aux2(a2DimX, a2LbX, a2UbX, a2EtaX, tau,
                a2EtaRatioX);

    vector<System*> auxs;
    auxs.push_back(&aux0);
    auxs.push_back(&aux1);
    auxs.push_back(&aux2);

    GBuchi abs("whirlpoolDiscrete.txt");

    TicToc tt;
    tt.tic();
    abs.initialize(&base, auxs);
    clog << "-------------------------------------------------------initialize: ";
    tt.toc();

//    abs.initializeProdGs(prod0AddG, 0);
//    abs.initializeProdGs(prod1AddG, 1);
//    abs.initializeProdGs(prod2AddG, 2);

//    abs.initializeBaseOs(baseAddO);

//    //TicToc tt;
//    tt.tic();
//    abs.computeBaseAbstractions(baseSysNext, baseRadNext, bx, bu);
//    abs.computeAuxAbstractions(aux0SysNext, aux0RadNext, a0x, au, 0);
//    abs.computeAuxAbstractions(aux1SysNext, aux1RadNext, a1x, au, 1);
//    abs.computeAuxAbstractions(aux2SysNext, aux2RadNext, a2x, au, 2);
//    abs.composeAbstractions();
//    clog << "-------------------------------------------------------abstraction process: ";
//    tt.toc();

//    int startAbs = 1-1;
//    int minToGoCoarser = 3;
//    int minToBeValid = 3;

//    abs.gBuchi(startAbs, minToGoCoarser, minToBeValid);
}


int main() {
    composition();
}
