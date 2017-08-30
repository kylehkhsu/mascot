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

// b means base, a means auxiliary
#define bDimX 2
#define bDimU 2
#define a0DimX 2
#define a1DimX 2
#define a2DimX 2
#define a3DimX 2
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

auto aux3SysNext = [](a3X_t &a3x, aU_t &au, double tau, OdeSolver solver) -> void {
    a3x[1] = a3x[1] + 1 * tau/0.6;
    if (a3x[1] > 2.75) {
        a3x[1] = -2.75;
    }
};

auto aux3RadNext= [](a3X_t &a3r, aU_t &au, double tau, OdeSolver solver) -> void {
    a3r[0] = 0;
    a3r[1] = 0;
};

// obstacle sets are only in terms of base
auto baseAddO = [](SymbolicSet* O) -> void {
    ;
};

// goal sets are only in terms of product systems
// ba0x = base-aux0 inflated state
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

auto prod3AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba3x)->bool {
        return sqrt(pow(ba3x[0]-ba3x[2], 2) + pow(ba3x[1]-ba3x[3], 2)) <= 1;
    };
    G->addByFunction(f);
};


void composition() {
    // 1. define all parameters for abstractions

    // base
    double bLbX[bDimX]={-6, -6};
    double bUbX[bDimX]={ 6,  6 };
    double bEtaX[bDimX]= {0.6, 0.6};
    double bLbU[bDimU]= {-2, 0.5};
    double bUbU[bDimU]= { 2,   1};
    double bEtaU[bDimU]= {0.5, 0.2};
    bX_t bx;
    bU_t bu;

    double bEtaRatioX[bDimX] = {2, 2};

    // common aux input

    aU_t au;

    // aux0
    double a0LbX[a0DimX] = {  -3, 0};
    double a0UbX[a0DimX] = {   3, 8};
    double a0EtaX[a0DimX] = {  1.0, 8};
    a0X_t a0x;
    double a0EtaRatioX[a0DimX] = {2, 1};

    // aux1
    double a1LbX[a1DimX] = {   0, -3};
    double a1UbX[a1DimX] = {   8,  3};
    double a1EtaX[a1DimX] = {  8,  1.0};
    a1X_t a1x;
    double a1EtaRatioX[a1DimX] = {1, 2};

    // aux2
    double a2LbX[a2DimX] = {  -3, -8};
    double a2UbX[a2DimX] = {   3,  0};
    double a2EtaX[a2DimX] = {  1.0,  8};
    a2X_t a2x;
    double a2EtaRatioX[a2DimX] = {2, 1};

    // aux3
    double a3LbX[a3DimX] = {  -8,  -3};
    double a3UbX[a3DimX] = {   0,   3};
    double a3EtaX[a3DimX] = {  8,   1.0};
    a3X_t a3x;
    double a3EtaRatioX[a3DimX] = {1, 2};

    // 2. Create System objects
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

    // 2. b) Create vector of System objects for auxs
    vector<System*> auxs;
    auxs.push_back(&aux0);
    auxs.push_back(&aux1);
    auxs.push_back(&aux2);
    auxs.push_back(&aux3);

    // 3. Construct GBuchi object with log file name
    GBuchi abs("whirlpool4P2A.txt");

    // 4. Initialize GBuchi object with base System and aux Systems vector
    abs.initialize(&base, auxs);

    // 5. Construct specification sets with (function pointer, auxiliary system index)
    abs.initializeProdGs(prod0AddG, 0);
    abs.initializeProdGs(prod1AddG, 1);
    abs.initializeProdGs(prod2AddG, 2);
    abs.initializeProdGs(prod3AddG, 3);

    // 6. Construct base system obstacle sets
    abs.initializeBaseOs(baseAddO);

    // 7. Compute base system abstractions
    abs.computeBaseAbstractions(baseSysNext, baseRadNext, bx, bu);

    // 8. Compute auxiliary systems abstractions
    abs.computeAuxAbstractions(aux0SysNext, aux0RadNext, a0x, au, 0);
    abs.computeAuxAbstractions(aux1SysNext, aux1RadNext, a1x, au, 1);
    abs.computeAuxAbstractions(aux2SysNext, aux2RadNext, a2x, au, 2);
    abs.computeAuxAbstractions(aux3SysNext, aux3RadNext, a3x, au, 3);

    // 8. b) Create composed/product transition systems
    abs.composeAbstractions();

    int startAbs = 1;
    int minToGoCoarser = 3;
    int minToBeValid = 3;

    // 9. Run synthesis algorithm
    abs.gBuchi(startAbs, minToGoCoarser, minToBeValid);
}


int main() {
    composition();
}
