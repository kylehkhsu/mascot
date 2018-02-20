#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsSafe.hh"

using namespace std;
using namespace scots;

#define dimX 3 /* dimension of the state space */
#define dimU 2 /* dimension of the input space */

/* bound on disturbance */
const double w[dimX] = {0.05, 0.05};

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* solution of the differential equations representing system dynamics (the result is stored in x)  */
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

/* the function which defines the safe set (either as a polytope or as an ellipsoid) */
auto unicycleAddS = [](SymbolicSet* S) -> void {
    double H[4*3]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};

    double h1[4] = {-2.1, 2.3, 0, 1.2};
    S->addGridPoints(); /* this action initializes the symbolic set S with the full state space */
    S->remPolytope(4, H, h1, OUTER); /* this action removes the polytope (the unsafe sets) defined by H and h from the state space */
    /* this is an example where we define the safe set by removing all the unsafe sets from the state space. Another option would be to *add* safe sets to the empty symbolic set S. The function to use in that case is S->addPolytope(4, H, h1, INNER). OUTER/INNER indicates the type of approximation used, which for safety should always be INNER for the safe sets, or OUTER for the unsafe sets. */
    ;
};

int main() {
    /* bounds on the state variables; 'lb', 'ub' stand for lower bound and upper boound resp. */
    double lbX[dimX] = {0, 0, -M_PI-0.4};
    double ubX[dimX] = {6.0, 1.8, M_PI+0.4};

    /* bounds on the input variables; 'lb', 'ub' stand for lower bound and upper boound resp. */
    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    
    /* discretization parameters of the *COARSEST* layer */
    double etaU[dimU] = {.3, .2};
    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    /* number of steps in numerical sollution of the ODEs */
    int nSubInt = 5;

    /* ratio of discretization parameters across adjacent layers (same for all adjacent layers)
     * the ratios must be multiple of 2 */
    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;

    /* total number of abstraction layers */
    int numAbs = 3;

    X_type x;
    U_type u;

    /* definition and initialization of different classes */
    System unicycle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);
    AdaptAbsSafe abs("unicycle_small_3A_safe.log"); /* the parameter is the name of the log file which contains many useful information */
    abs.initialize(&unicycle, unicycleAddS);

    TicToc tt_tot;
    tt_tot.tic();
    abs.onTheFlySafe(sysNext, radNext, x, u); /* lazy abstraction + synthesis happen here */
    clog << "------------------------------------Total time:";
    tt_tot.toc();
    
    /* DISCLAIMER: there is an unresolved bug which causes an exception in the end. For the moment, pleasse ignore the message, as this doesn't cause any difference to the results. */
}



