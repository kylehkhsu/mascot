#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "cuddObj.hh"

#include "Adaptive.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

using namespace std;
using namespace scots;


/* state space dim */
int dimX = 2;
int dimU = 2;

/* data types for the ode solver */
typedef std::array<double,2> X_type;
typedef std::array<double,2> U_type;

/* sampling time */
double tau = 1;
/* number of intermediate steps in the ode solver */
int nint = 5;
OdeSolver ode_solver(dimX,nint,tau);

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u) -> void {

//    /* the ode describing the simple */
//    auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
//        xx[0] = u[0];
//    };
//    ode_solver(rhs,x,u);
    x[0] += u[0];
    x[1] += u[1];

};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u) -> void {
    r[0] = 0.6;
    r[1] = 0.6;
};

// forward declarations
SymbolicSet simpleCreateStateSpace(Cudd &mgr);
SymbolicSet simpleCreateInputSpace(Cudd &mgr);

int main() {
    double lbX[dimX]={-10, -10};
    double ubX[dimX]={10, 10};
    double etaX[dimX]={0.5, 0.5};
    double lbU[dimU]={-1, -1};
    double ubU[dimU]={1, 1};
    double etaU[dimU]={1, 1};

    double q[dimX] = {7, 8};
    double H[dimX * dimX]={0.65, 0,
                 0, 0.65};
    double c[dimX] = {-7, -7};
    Cudd mgr;

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX,
                                 dimU, lbU, ubU, etaU,
                                 q, H, c);
    abs.iteration(sysNext, radNext);
//    abs.test();

//    Cudd mgr;
//    TicToc tt;

//    SymbolicSet X = simpleCreateStateSpace(mgr);
//    X.writeToFile("X.bdd");
//    X.complement();
//    X.writeToFile("X_obst.bdd");
//    X.complement();

//    SymbolicSet U = simpleCreateInputSpace(mgr);
//    cout << endl << "Input space details:" << endl;
//    U.printInfo(1);

//    SymbolicSet G (X);
//    double H[9]={0.65, 0,
//                  0, 0.65};
//    double c[3] = {-7, -7};
//    G.addEllipsoid(H,c, INNER);
////    G.addPoint(g0);

//    G.setSymbolicSet(G.getSymbolicSet() & X.getSymbolicSet());
//    G.writeToFile("G.bdd");

//    SymbolicSet X2(X,1);
//    SymbolicModelGrowthBound<X_type,U_type> abstraction(&X, &U, &X2);
//    tt.tic();
//    abstraction.computeTransitionRelation(sysNext, radNext);
//    cout << endl;
//    tt.toc();
//    cout << endl << "Number of elements in the transition relation: " << abstraction.getSize() << endl;

//    int verbose = 1;
//    FixedPoint fp(&abstraction);
//    BDD T = G.getSymbolicSet();
//    tt.tic();
//    BDD C = fp.reach(T,verbose);
//    tt.toc();
//    SymbolicSet controller(X,U);
//    controller.setSymbolicSet(C);
//    controller.writeToFile("C.bdd");

//    SymbolicSet tr = abstraction.getTransitionRelation();
//    tr.writeToFile("T.bdd");
//    return 1;
}

SymbolicSet simpleCreateStateSpace(Cudd &mgr) {

    /* setup the workspace of the synthesis problem and the uniform grid */
    /* lower bounds of the hyper rectangle */
    double lb[dimX]={-10, -10};
    /* upper bounds of the hyper rectangle */
    double ub[dimX]={10, 10};
    /* grid node distance diameter */
    double eta[dimX]={0.5, 0.5};

    /* eta is added to the bound so as to ensure that the whole
   * [0,10]x[0,10]x[-pi-eta,pi+eta] is covered by the cells */

    SymbolicSet ss(mgr,dimX,lb,ub,eta);

    /* add the grid points to the SymbolicSet ss */
    ss.addGridPoints();

    double H[4*dimX]={-1, 0,
                      1, 0,
                      0,-1,
                      0, 1};
    double h1[4] = {5, 4, 2, 5};
    ss.remPolytope(4, H, h1, OUTER);
    double h2[4] = {10, -6, 2, 5};
    ss.remPolytope(4, H, h2, OUTER);
    double h3[4] = {2, 8, 7, -4};
    ss.remPolytope(4, H, h3, OUTER);


    return ss;
}

SymbolicSet simpleCreateInputSpace(Cudd &mgr) {

    /* lower bounds of the hyper rectangle */
    double lb[dimU]={-1, -1};
    /* upper bounds of the hyper rectangle */
    double ub[dimU]={1, 1};
    /* grid node distance diameter */
    double eta[dimU]={1, 1};

    SymbolicSet is(mgr,dimU,lb,ub,eta);
    is.addGridPoints();

    return is;
}
