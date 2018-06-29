#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsReach.hh"
// #include "Reach.hh"

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
auto radNext = [](X_type &x, X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
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

// auto simpleAddG = [](SymbolicSet* G) -> void {
//     double H[4*2]={-1, 0,
//                     1, 0,
//                     0,-1,
//                     0, 1};
//     double h1[4] = {-2.5, 3.5, 1, 0.7};
//     G->addPolytope(4, H, h1, INNER);
// };

auto target = [](const scots::abs_type &abs_state, const scots::UniformGrid* ss) {
  X_type x;
  ss->itox(abs_state,x);
  std::vector<double> s_eta = ss->get_eta();
  /* function returns 1 if cell associated with x is in target set  */
  if (2.5 <= (x[0]) && (x[0]) <= 3.5 &&
      -0.7 <= (x[1]) && (x[1]) <= 0.7)
    return true;
  return false;
};

// auto simpleAddO = [](SymbolicSet* O) -> void {
//     double H[4*2]={-1, 0,
//                     1, 0,
//                     0,-1,
//                     0, 1};
//
//     double h1[4] = {-2.1, 2.3, 0, 1.2};
//     O->addPolytope(4, H, h1, OUTER);
//     ;
// };

auto obstacle = [](const scots::abs_type &abs_state, const scots::UniformGrid &ss) {
  X_type x;
  ss.itox(abs_state,x);
  if (2.1 <= (x[0]) && (x[0]) <= 2.3 &&
      0 <= (x[1]) && (x[1]) <= 1.2)
    return true;
  return false;
};

// auto simpleAddI = [](SymbolicSet* I) -> void {
//     double q[3] = {0.25, 0.25, 0};
//     I->addPoint(q);
// };

int main() {

    double lbX[dimX]={0, 0};
    double ubX[dimX]={6, 1.8};

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
    int numAbs = 1;
    bool lazy = true;
    int readAbs = 0;
    int p_ = 2;
    bool verbose = true;

    System system(dimX, lbX, ubX, etaX, tau,
                  dimU, lbU, ubU, etaU,
                  etaRatio, tauRatio, nSubInt, numAbs);

    AdaptAbsReach syn("simple_small3A.log");
    syn.initialize(&system, target, obstacle, verbose);

    TicToc tt_tot;
    tt_tot.tic();
    syn.onTheFlyReach(p_, sysNext, radNext, x, u, lazy, readAbs);
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
