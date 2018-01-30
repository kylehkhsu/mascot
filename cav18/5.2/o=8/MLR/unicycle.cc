#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "UpfrontReach.hh"

using namespace std;
using namespace scots;
using namespace helper;

#define dimX 3
#define dimU 2

const double w[dimX] = {0.05, 0.05};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* we integrate the unicycle ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](X_type& xx,  const X_type &x, U_type &u) -> void {
      xx[0] = u[0]*std::cos(x[2]);
      xx[1] = u[0]*std::sin(x[2]);
      xx[2] = u[1];
  };
  solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */

auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0] + (r[2]*std::abs(u[0]) + w[0]) * tau;
    r[1] = r[1] + (r[2]*std::abs(u[0]) + w[1]) * tau;
};

auto unicycleAddG = [](SymbolicSet* G) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
    double h1[4] = {-22.1, 22.9, -5.3, 6.1};
    G->addPolytope(4, H, h1, INNER);
};

auto unicycleAddO = [](SymbolicSet* O) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double wall1[4] = {-1.9, 2.3, -1.9, 6};
    O->addPolytope(4, H, wall1, OUTER);

    double wall2[4] = {-4.3, 4.7, -0, 4.1};
    O->addPolytope(4, H, wall2, OUTER);

    double wall3[4] = {-6.7, 7.1, -1.9, 6};
    O->addPolytope(4, H, wall3, OUTER);

    double wall4[4] = {-9.1, 9.5, -0, 4.1};
    O->addPolytope(4, H, wall4, OUTER);

    double wall5[4] = {-11.5, 11.9, -1.9, 6};
    O->addPolytope(4, H, wall5, OUTER);

    double wall6[4] = {-13.9, 14.3, -0, 4.1};
    O->addPolytope(4, H, wall6, OUTER);

    double wall7[4] = {-16.3, 16.7, -1.9, 6};
    O->addPolytope(4, H, wall7, OUTER);

    double wall8[4] = {-18.7, 19.1, -0, 4.1};
    O->addPolytope(4, H, wall8, OUTER);

    double wall9[4] = {-21.1, 21.5, -1.9, 6};
    O->addPolytope(4, H, wall9, OUTER);

    double obs1[4] = {-0.7, 1.1, -3.7, 4.1};
    O->addPolytope(4, H, obs1, OUTER);

    double obs2[4] = {-0.7, 1.1, -1.9, 2.3};
    O->addPolytope(4, H, obs2, OUTER);

    double obs3[4] = {-3.1, 3.5, -1.9, 2.3};
    O->addPolytope(4, H, obs3, OUTER);

    double obs4[4] = {-3.1, 3.5, -3.7, 4.1};
    O->addPolytope(4, H, obs4, OUTER);

    double obs5[4] = {-5.5, 5.9, -3.7, 4.1};
    O->addPolytope(4, H, obs5, OUTER);

    double obs6[4] = {-5.5, 5.9, -1.9, 2.3};
    O->addPolytope(4, H, obs6, OUTER);

    double obs7[4] = {-7.9, 8.3, -1.9, 2.3};
    O->addPolytope(4, H, obs7, OUTER);

    double obs8[4] = {-7.9, 8.3, -3.7, 4.1};
    O->addPolytope(4, H, obs8, OUTER);
};

auto unicycleAddI = [](SymbolicSet* I) -> void {
    double q[3] = {0.55, 5.85, -M_PI/2};
    I->addPoint(q);
};

int main() {

    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {23.4, 6, M_PI+0.4};

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    int nSubInt = 5;

    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;

    int numAbs = 3;

    X_type x;
    U_type u;

    System unicycle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    int m = 2;
    UpfrontReach abs("unicycle.log");
    abs.initialize(&unicycle, 0, unicycleAddO);
    abs.initializeReach(unicycleAddG, unicycleAddI);

    TicToc timer;
    timer.tic();
    abs.computeAbstractions(sysNext, radNext, x, u);
    abs.upfrontReach(m);
    clog << "------------------------------------Total time: " << timer.toc() << " seconds.\n";
}



