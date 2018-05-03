#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "Adaptive.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define dimX 2
#define dimU 2

auto sandboxAddO = [](SymbolicSet* O) -> void {
    ;
};

void testFiner() {
    double lbX[dimX] = {-3, -3};
    double ubX[dimX] = { 3,  3};
    double etaX[dimX] = {1, 1};
    double tau = 1;

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    double etaRatio[dimX] = {2, 1};
    double tauRatio = 2;

    int nSubInt = 5;
    int numAbs = 3;

    int readAbs = 0;

    System sandbox(dimX, lbX, ubX, etaX, tau,
                   dimU, lbU, ubU, etaU,
                   etaRatio, tauRatio, nSubInt, numAbs);
    Adaptive abs("sandboxFiner.txt");
    abs.initialize(&sandbox, readAbs, sandboxAddO);

    SymbolicSet Zc(*abs.Xs_[1]);

    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {0, 2, 1, 1};

    Zc.addPolytope(4, H, c, OUTER);
    SymbolicSet Zf(*abs.Xs_[2]);

    abs.finer(&Zc, &Zf, 1);

    abs.Xs_[0]->writeToFile("X.bdd");
    Zf.writeToFile("Zf.bdd");
    Zc.writeToFile("Zc.bdd");
}

void testCoarser() {
    double lbX[dimX] = {-2, -2};
    double ubX[dimX] = { 2,  2};
    double etaX[dimX] = {1, 1};
    double tau = 1;

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    double etaRatio[dimX] = {2, 1};
    double tauRatio = 2;

    int nSubInt = 5;
    int numAbs = 3;

    int readAbs = 0;

    System sandbox(dimX, lbX, ubX, etaX, tau,
                   dimU, lbU, ubU, etaU,
                   etaRatio, tauRatio, nSubInt, numAbs);
    Adaptive abs("sandboxCoarser.txt");
    abs.initialize(&sandbox, readAbs, sandboxAddO);

    SymbolicSet Zf(*abs.Xs_[2]);

    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {0, 2, 1, 1};

    Zf.addPolytope(4, H, c, OUTER);
    SymbolicSet Zc(*abs.Xs_[1]);

    abs.coarser(&Zc, &Zf, 1);

    abs.Xs_[0]->writeToFile("X.bdd");
    Zf.writeToFile("Zf.bdd");
    Zc.writeToFile("Zc.bdd");

    abs.Xs_[0]->addGridPoints();
    abs.Xs_[0]->printInfo(2);
}

void hardCode() {
    Cudd mgr;
    double lb[dimX-1] = {-3.01};
    double ub[dimX-1] = { 3.01};
    double etaf[dimX-1] = {2};

    SymbolicSet Xf(mgr, dimX-1, lb, ub, etaf, 0);
    Xf.addGridPoints();
//    Xf.symbolicSet_ = mgr.bddOne();


    double q[dimX-1] = {-1};
//    Xf.remPoint(q);
//    q[0] = 1;
//    Xf.addPoint(q);
//    q[1] = -1;
//    Xf.addPoint(q);
//    q[0] = -1;
//    Xf.addPoint(q);

    double gridPoints[4] = {0};
    Xf.copyGridPoints(gridPoints);
    printArray(gridPoints, 4);

    Xf.printInfo(2);


//    SymbolicSet Zf(Xf);

//    double H[4*2]={-1, 0,
//                      1, 0,
//                      0,-1,
//                      0, 1};
//    double c[4] = {0, 2, 1, 1};
//    Zf.addPolytope(4, H, c, OUTER);

//    Zf.printInfo(2);

    Xf.writeToFile("Xf.bdd");

    SymbolicSet Xread(mgr, "Xf.bdd");
    Xread.printInfo(2);
//    Zf.writeToFile("Zf.bdd");


//    SymbolicSet Tf(Xf);
//    double q[dimX] = {-2.25, -2.25};
//    Tf.addPoint(q);

//    Tf.printInfo(2);


//    BDD* vars = new BDD[2];
//    vars[0] = mgr.bddVar(Xf.idBddVars_[0]);
//    vars[1] = mgr.bddVar(Xf.idBddVars_[4]);

//    BDD cube = mgr.bddComputeCube(vars, NULL, 2);
//    SymbolicSet cubeSS(Xf);
//    cubeSS.symbolicSet_ = cube;
//    cout << "cube:\n";
//    cubeSS.printInfo(2);

//    SymbolicSet Zf2(Zf);
////    Zf2.symbolicSet_ = !(!Zf.symbolicSet_.ExistAbstract(cube));
//    Zf2.symbolicSet_ = Zf.symbolicSet_.UnivAbstract(cube);
//    Zf2.printInfo(2);
//    Zf2.writeToFile("Zf2.bdd");

//    delete[] vars;



//    double etac[dimX] = {1, 1};

//    SymbolicSet Xc(mgr, dimX, lb, ub, etac, 0);
//    Xc.addGridPoints();

//    int permute[Xf.nvars_] = {0, 8, 9, 10, 4, 11, 12, 13};

//    SymbolicSet Zc(Xc);
//    Zc.symbolicSet_ = Zf2.symbolicSet_.Permute(permute);

//    Zc.printInfo(1);

//    Zc.writeToFile("Zc.bdd");
}

int main() {
//    testCoarser();
//    testFiner();
    hardCode();

}