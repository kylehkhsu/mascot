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


int main() {

    Cudd mgr;
    double lb[dimX] = {-3, -3};
    double ub[dimX] = { 3,  3};
    double etaf[dimX] = {0.5, 0.5};

    SymbolicSet Xf(mgr, dimX, lb, ub, etaf, 0);
    Xf.addGridPoints();

    SymbolicSet Zf(Xf);

    double H[4*2]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double c[4] = {0, 2, 1, 1};
    Zf.addPolytope(4, H, c, OUTER);

    Zf.printInfo(2);

    Xf.writeToFile("Xf.bdd");
    Zf.writeToFile("Zf.bdd");


    SymbolicSet Tf(Xf);
    double q[dimX] = {-2.25, -2.25};
    Tf.addPoint(q);

    Tf.printInfo(2);


    BDD* vars = new BDD[2];
    vars[0] = mgr.bddVar(Xf.idBddVars_[0]);
    vars[1] = mgr.bddVar(Xf.idBddVars_[4]);

    BDD cube = mgr.bddComputeCube(vars, NULL, 2);
    SymbolicSet cubeSS(Xf);
    cubeSS.symbolicSet_ = cube;
    cout << "cube:\n";
    cubeSS.printInfo(2);

    SymbolicSet Zf2(Zf);
//    Zf2.symbolicSet_ = !(!Zf.symbolicSet_.ExistAbstract(cube));
    Zf2.symbolicSet_ = Zf.symbolicSet_.UnivAbstract(cube);
    Zf2.printInfo(2);
    Zf2.writeToFile("Zf2.bdd");

    delete[] vars;



    double etac[dimX] = {1, 1};

    SymbolicSet Xc(mgr, dimX, lb, ub, etac, 0);
    Xc.addGridPoints();

    int permute[Xf.nvars_] = {0, 8, 9, 10, 0, 11, 12, 13};

    SymbolicSet Zc(Xc);
    Zc.symbolicSet_ = Zf2.symbolicSet_.Permute(permute);

    Zc.printInfo(1);

    Zc.writeToFile("Zc.bdd");

}
