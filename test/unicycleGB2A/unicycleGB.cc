// unicycleGB with stationary auxiliary targets
#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "GBuchi.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define bDimX 3
#define bDimU 2
#define a0DimX 2
#define a1DimX 2
#define a2DimX 2
#define aDimU 1

const int numAbs = 3-1;
const double tau = 0.9/2;
const double tauRatio = 2;
const int nSubInt = 5;

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
        dbxdt[0] = bu[0]*std::cos(bx[2]);
        dbxdt[1] = bu[0]*std::sin(bx[2]);
        dbxdt[2] = bu[1];
    };
    solver(ODE, bx, bu);
};

auto baseRadNext = [](bX_t &br, bU_t &bu, double tau, OdeSolver solver) -> void {
    br[0] = br[0] + (br[2]*std::abs(bu[0]) + bw[0]) * tau;
    br[1] = br[1] + (br[2]*std::abs(bu[0]) + bw[1]) * tau;
};

auto aux0SysNext = [](a0X_t &a0x, aU_t &au, double tau, OdeSolver solver) -> void {
    ;
};

auto aux0RadNext= [](a0X_t &a0r, aU_t &au, double tau, OdeSolver solver) -> void {
    a0r[0] = 0;
    a0r[1] = 0;
};

auto aux1SysNext = [](a1X_t &a1x, aU_t &au, double tau, OdeSolver solver) -> void {
    ;
};

auto aux1RadNext= [](a1X_t &a1r, aU_t &au, double tau, OdeSolver solver) -> void {
    a1r[0] = 0;
    a1r[1] = 0;
};

auto aux2SysNext = [](a2X_t &a2x, aU_t &au, double tau, OdeSolver solver) -> void {
    ;
};

auto aux2RadNext= [](a2X_t &a2r, aU_t &au, double tau, OdeSolver solver) -> void {
    a2r[0] = 0;
    a2r[1] = 0;
};

auto baseAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double h1[4] = {-1.9, 2.3, -1.9, 12};
    O->addPolytope(4, H, h1, OUTER);

    double h2[4] = {-4.3, 4.7, -0, 10.1};
    O->addPolytope(4, H, h2, OUTER);

    double h3[4] = {-6.7, 7.1, -1.9, 12};
    O->addPolytope(4, H, h3, OUTER);

    double h6[4] = {-9.1, 9.5, -0, 10.1};
    O->addPolytope(4, H, h6, OUTER);

    double h4[4] = {-2.5, 3.2, -3.7, 4.6};
    O->addPolytope(4, H, h4, OUTER);

    double h5[4] = {-5.39, 6.5, -4.9, 6.5};
    O->addPolytope(4, H, h5, OUTER);
    ;
};

auto prod0AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba0x)->bool {
        return abs(ba0x[0]-ba0x[3]) <= 0.2 && abs(ba0x[1]-ba0x[4]) <= 0.2;
    };
    G->addByFunction(f);
};

auto prod1AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba1x)->bool {
        return abs(ba1x[0]-ba1x[3]) <= 0.2 && abs(ba1x[1]-ba1x[4]) <= 0.2;
    };
    G->addByFunction(f);
};

auto prod2AddG = [](SymbolicSet* G) -> void {
    auto f = [](double* ba2x)->bool {
        return abs(ba2x[0]-ba2x[3]) <= 0.2 && abs(ba2x[1]-ba2x[4]) <= 0.2;
    };
    G->addByFunction(f);
};


void composition() {
    double bLbX[bDimX]={0, 0, -M_PI-0.4};
    double bUbX[bDimX]={11.4, 11.4, M_PI+0.4};
    double bEtaX[bDimX]= {0.6/2, 0.6/2, 0.3/2};
    double bLbU[bDimU]= {-1.2, -1.6};
    double bUbU[bDimU]= {1.2, 1.6};
    double bEtaU[bDimU]= {.3, .2};
    bX_t bx;
    bU_t bu;

    double bEtaRatioX[bDimX] = {2, 2, 2};

    aU_t au;

    double a0LbX[a0DimX] = {0, 0};
    double a0UbX[a0DimX] = {0.7, 22};
    double a0EtaX[a0DimX] = { 0.7, 22};
    a0X_t a0x;
    double a0EtaRatioX[a0DimX] = {1, 1};

    double a1LbX[a1DimX] = {0, 0};
    double a1UbX[a1DimX] = {22, 0.7};
    double a1EtaX[a1DimX] = { 22, 0.7};
    a1X_t a1x;
    double a1EtaRatioX[a1DimX] = {1, 1};

    double a2LbX[a2DimX] = {0, 0};
    double a2UbX[a2DimX] = {12.4, 22};
    double a2EtaX[a2DimX] = {  12.4, 22};
    a2X_t a2x;
    double a2EtaRatioX[a2DimX] = {1, 1};
    
    int startAbs = 1-1;
//    int readAbs = 0;

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

    GBuchi abs("unicycleGB2A.txt");

    TicToc tt;
    tt.tic();
    abs.initialize(&base, auxs);
    clog << "-------------------------------------------------------initialize: ";
    tt.toc();

    abs.initializeProdGs(prod0AddG, 0);
    abs.initializeProdGs(prod1AddG, 1);
    abs.initializeProdGs(prod2AddG, 2);

    abs.initializeBaseOs(baseAddO);

    abs.verify();

    tt.tic();
    abs.computeBaseAbstractions(baseSysNext, baseRadNext, bx, bu);
    abs.computeAuxAbstractions(aux0SysNext, aux0RadNext, a0x, au, 0);
    abs.computeAuxAbstractions(aux1SysNext, aux1RadNext, a1x, au, 1);
    abs.computeAuxAbstractions(aux2SysNext, aux2RadNext, a2x, au, 2);
    abs.composeAbstractions();
    clog << "-------------------------------------------------------abstraction process: ";
    tt.toc();

    int minToGoCoarser = 5;
    int minToBeValid = 5;

    abs.gBuchi(startAbs, minToGoCoarser, minToBeValid);
}


int main() {
    composition();
}
