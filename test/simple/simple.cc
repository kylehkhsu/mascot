#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsReach.hh"
#include "Reach.hh"

using namespace std;
using namespace scots;

/* state space dim */
#define dimX 2
#define dimU 2

/* disturbance */
const double w[dimX] = {0.05, 0.05};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    x[0] += u[0] * tau;
    x[1] += u[1] * tau;
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        double L[2][2];
        L[0][0] = 0;
        L[0][1] = 0;
        L[1][0] = 0;
        L[1][1] = 0;

        drdt[0] = L[0][0]*r[0] + L[0][1]*r[1] + w[0];
        drdt[1] = L[1][0]*r[0] + L[1][0]*r[1] + w[1];
    };
    solver(radODE, r, u);
};

auto simpleAddG = [](SymbolicSet* G) -> void {
    double H[4*2]={-1, 0,
                    1, 0,
                    0,-1,
                    0, 1};
    double h1[4] = {-9.5, 11.4, 0, 3};
    G->addPolytope(4, H, h1, INNER);
};

auto simpleAddO = [](SymbolicSet* O) -> void {
    double H[4*2]={-1, 0,
                    1, 0,
                    0,-1,
                    0, 1};

    double h1[4] = {-1.9, 2.3, -1.9, 12};
    O->addPolytope(4, H, h1, OUTER);

    double h2[4] = {-4.3, 4.7, -0, 10.1};
    O->addPolytope(4, H, h2, OUTER);

    double h3[4] = {-6.7, 7.1, -1.9, 12};
    O->addPolytope(4, H, h3, OUTER);

    double h4[4] = {-2.5, 3.2, -3.7, 4.6};
    O->addPolytope(4, H, h4, OUTER);

    double h5[4] = {-5.39, 6.5, -4.9, 6.5};
    O->addPolytope(4, H, h5, OUTER);

    double h6[4] = {-7, 8.5, -7, 9};
    O->addPolytope(4, H, h6, OUTER);

    double h7[4] = {-0.5, 2, -6, 8};
    O->addPolytope(4, H, h7, OUTER);

    double h8[4] = {-9.1, 9.5, -0, 10.4};
    O->addPolytope(4, H, h8, OUTER);
};

auto simpleAddI = [](SymbolicSet* I) -> void {
    double q[3] = {0.55, 10.85, -M_PI/2};
    I->addPoint(q);
};

int main() {

    double lbX[dimX]={0, 0};
    double ubX[dimX]={11.4, 11.4};

    double lbU[dimU]={-1.3, -1.3};
    double ubU[dimU]= {1.3, 1.3};
    double etaU[dimU]= {0.5, 0.5};

    X_type x;
    U_type u;

    double etaRatio[dimX] = {2, 2};
    double tauRatio = 2;
    int nSubInt = 5;

    double etaX[dimX]= {0.6, 0.6};
    double tau = 0.9;
    int numAbs = 3;

    System system(dimX, lbX, ubX, etaX, tau,
                  dimU, lbU, ubU, etaU,
                  etaRatio, tauRatio, nSubInt, numAbs);

    AdaptAbsReach syn("simple3A.log");
    syn.initialize(&system, simpleAddO, simpleAddG);

    TicToc tt_tot;
    tt_tot.tic();
    syn.onTheFlyReach(sysNext, radNext, x, u);
    clog << "------------------------------------Total time:";
    tt_tot.toc();

//    int startAbs = 0;
//    int minToGoCoarser = 1;
//    int minToBeValid = 2;
//    int earlyBreak = 0;

//    Reach abs("simple3Aprev.log");
//    abs.initialize(&system, 0, simpleAddO);
//    abs.initializeReach(simpleAddG, simpleAddI);

//    TicToc tt_tot;
//    tt_tot.tic();
//    abs.computeAbstractions(sysNext, radNext, x, u);
//    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
//    clog << "------------------------------------Total time:";
//    tt_tot.toc();
}
