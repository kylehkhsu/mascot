#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsReach.hh"
//#include "Reach.hh"

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

auto radNext = [](X_type &x, X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    r[0] = r[0] + (r[2]*std::abs(u[0]) + w[0]) * tau;
    r[1] = r[1] + (r[2]*std::abs(u[0]) + w[1]) * tau;
};

auto target = [](const scots::abs_type &abs_state, const scots::UniformGrid* ss) {
	// X_type t_lb = { { 4.5, 0.1, -M_PI - 0.8 } }; // 0.8 offset is made up
	// X_type t_ub = { { 6.0, 0.61, M_PI + 0.8 } }; // 0.8 offset is made up
  X_type t_lb = { { 4.5, 0.1, -3.3 } }; // 0.8 offset is made up
	X_type t_ub = { { 6.0, 0.61, 3.3 } }; // 0.8 offset is made up
	X_type c_lb;
	X_type c_ub;
	X_type z = { {0, 0, 0} }; // assume no measurement error
	/* center of cell associated with abs_state is stored in x */
	X_type x;
	ss->itox(abs_state, x);
	/* hyper-interval of the quantizer symbol with perturbation */
	int dim = ss->get_dim();
	std::vector<double> etaX = ss->get_eta();
	for (int i = 0; i<dim; i++) {
		c_lb[i] = x[i] - etaX[i] / 2.0 - z[i];
		c_ub[i] = x[i] + etaX[i] / 2.0 + z[i];
	}
	if (t_lb[0] <= c_lb[0] && c_ub[0] <= t_ub[0] &&
		t_lb[1] <= c_lb[1] && c_ub[1] <= t_ub[1] &&
		t_lb[2] <= c_lb[2] && c_ub[2] <= t_ub[2]) {
		return true;
	}
	return false;

};

auto obstacle = [](const scots::abs_type &abs_state, const scots::UniformGrid &ss) {
	// X_type t_lb = { { 1.9, 0, -M_PI - 0.4 } };
	// X_type t_ub = { { 2.3, 1.2, M_PI + 0.4 } };
  X_type t_lb = { { 1.9, 0, -3.3 } };
	X_type t_ub = { { 2.3, 1.2, 3.3 } };
	X_type c_lb;
	X_type c_ub;
	X_type z = { { 0, 0, 0 } }; // assume no measurement error
	/* center of cell associated with abs_state is stored in x */
	X_type x;
	ss.itox(abs_state, x);
	/* hyper-interval of the quantizer symbol with perturbation */
	int dim = ss.get_dim();
	std::vector<double> etaX = ss.get_eta();
	for (int i = 0; i<dim; i++) {
		c_lb[i] = x[i] - etaX[i] / 2.0 - z[i];
		c_ub[i] = x[i] + etaX[i] / 2.0 + z[i];
	}
	if (t_lb[0] <= c_lb[0] && c_ub[0] <= t_ub[0] &&
		t_lb[1] <= c_lb[1] && c_ub[1] <= t_ub[1]) { //&&
		//t_lb[2] <= c_lb[2] && c_ub[2] <= t_ub[2]) {
		return true;
	}
	return false;

};

//auto unicycleAddI = [](SymbolicSet* I) -> void {
//    double q[3] = {0.5, 0.5, 0};
//    I->addPoint(q);
//};

int main() {

    // double lbX[dimX] = {0, 0, -M_PI-0.4};
    // double ubX[dimX] = {6.0, 1.8, M_PI+0.4};
    double lbX[dimX] = {0, 0, -3.3};
    double ubX[dimX] = {6.0, 1.8, 3.3};
    // double lbX[dimX] = {0, 0, -0.3};
    // double ubX[dimX] = {0.6, 0.6, 0.3};

    double lbU[dimU] = {-1.2, -1.6};
    double ubU[dimU] = {1.2, 1.6};
    double etaU[dimU] = {.3, .2};

    int nSubInt = 5;

    double etaX[dimX] = {0.6, 0.6, 0.3};
    double tau = 0.9;

    double etaRatio[dimX] = {2, 2, 2};
    double tauRatio = 2;

    int numAbs = 3;
    int readAbs = 0;
	int p_ = 2;
    bool lazy = false;

    X_type x;
    U_type u;

    System unicycle(dimX, lbX, ubX, etaX, tau,
                    dimU, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);

    AdaptAbsReach abs("unicycle3A.log");
    abs.initialize(&unicycle, target, obstacle);

    TicToc tt_tot;
    tt_tot.tic();
    abs.onTheFlyReach(p_, sysNext, radNext, x, u, lazy);
    clog << "------------------------------------Total time:";
    tt_tot.toc();

//    int startAbs = 0;
//    int minToGoCoarser = 1;
//    int minToBeValid = 2;
//    int earlyBreak = 0;
//
//    Reach abs("unicycle_small.log");
//    abs.initialize(&unicycle, readAbs, unicycleAddO);
//    abs.initializeReach(unicycleAddG, unicycleAddI);
//    abs.computeAbstractions(sysNext, radNext, x, u);
//
//    abs.reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);
}
