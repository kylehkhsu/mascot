/*
 * aircraft.cc
 *
 *  created: Oct 2016
 *  author: Matthias Rungger
 *  modified by Kaushik Mallik
 *
 */

/*
 * information about this example is given in
 * http://arxiv.org/abs/1503.03715
 * doi: 10.1109/TAC.2016.2593947
 * doi: 10.1109/CDC.2015.7403185
 *
 */

#include <iostream>
#include <array>
#define _USE_MATH_DEFINES

/* SCOTS header */
//#include "scots.hh"
/* ode solver */
//#include "RungeKutta4.hh"
#include "AdaptAbsReach.hh"

using namespace std;
using namespace scots;

/* time profiling */
//#include "TicToc.hh"
/* memory profiling */
//#include <sys/time.h>
//#include <sys/resource.h>
// struct rusage usage;

#define state_dim 3
#define input_dim 2
/* state space dim */
//const int state_dim=3;
///* input space dim */
//const int input_dim=2;
/* sampling time */
//const double tau = 0.25;

/* data types of the state space elements and input 
 * space elements used in uniform grid and ode solver */
//using state_type = std::array<double,state_dim>;
//using input_type = std::array<double,input_dim>;

/* data types for the ode solver */
typedef std::array<double, state_dim> state_type;
typedef std::array<double, input_dim> input_type;

/* we integrate the aircraft ode by 0.25 sec (the result is stored in x)  */
auto aircraft_post = [] (state_type &x, input_type &u, double tau, OdeSolver solver) -> void {
  /* the ode describing the aircraft */
  auto rhs =[] (state_type& xx,  const state_type &x, input_type &u) -> void {
    double mg = 60000.0*9.81;
    double mi = 1.0/60000;
    double c=(1.25+4.2*u[1]);
    xx[0] = mi*((u[0]-16000)*std::cos(u[1])-(2.7+3.08*c*c)*x[0]*x[0]-mg*std::sin(x[1])); // shift u[0] by 16000 to account for grid alignment offset
    xx[1] = (1.0/(60000*x[0]))*((u[0]-16000)*std::sin(u[1])+68.6*c*x[0]*x[0]-mg*std::cos(x[1]));
    xx[2] = x[0]*std::sin(x[1]);
  };
  /* use 10 intermediate steps */
  solver(rhs,x,u);
};

/* we integrate the growth bound by 0.25 sec (the result is stored in r)  */
auto radius_post = [] (state_type &r, input_type &u, double tau, OdeSolver solver) -> void {
  /* the ode for the growth bound */
  auto rhs =[] (state_type& rr,  const state_type &r, input_type &u) -> void {
    /* lipschitz matrix computed with mupad/mathematica check the ./helper directory */
    double L[3][2];
    L[0][0]=-0.00191867*(2.7+3.08*(1.25+4.2*u[1])*(1.25+4.2*u[1]));
    L[0][1]=9.81;
    L[1][0]=0.002933+0.004802*u[1];
    L[1][1]=0.003623;
    L[2][0]=0.07483;
    L[2][1]=83.22;
    /* to account for input disturbances */
    const state_type w={{.108,0.002,0}};
    rr[0] = L[0][0]*r[0]+L[0][1]*r[1]+w[0]; /* L[0][2]=0 */
    rr[1] = L[1][0]*r[0]+L[1][1]*r[1]+w[1]; /* L[1][2]=0 */
    rr[2] = L[2][0]*r[0]+L[2][1]*r[1]+w[2]; /* L[2][2]=0 */
  };
  /* use 10 intermediate steps */
  solver(rhs,r,u);
};

auto aircraftAddG = [](SymbolicSet* G) -> void {
    double H[6*3]={-1, 0, 0,
                    1, 0, 0,
                    0, -1, 0,
                    0, 1, 0,
                    0, 0, -1,
                    0, 0, 1};
    double h[6] = {-63, 75, 3*M_PI/180, 0, 0, 2.5};
    G->addPolytope(6, H, h, INNER);
};

auto aircraftAddO = [](SymbolicSet* O) -> void {
    // double H[6*3]={-1, 0, 0,
    //                 1, 0, 0,
    //                 0, -1, 0,
    //                 0, 1, 0,
    //                 0, 0, -1,
    //                 0, 0, 1};
    // double h[6] = {-58, 83, 3*M_PI/180, 0, 0, 56};
    // O->addGridPoints();
    // O->remPolytope(6, H, h, INNER);
    ;
};

auto aircraftAddI = [](SymbolicSet* I) -> void {
    // double H[6*3]={-1, 0, 0,
    //                 1, 0, 0,
    //                 0, -1, 0,
    //                 0, 1, 0,
    //                 0, 0, -1,
    //                 0, 0, 1};
    // double h[6] = {-80, 82, 2*M_PI/180, -1*M_PI/180, -54.9, 55.1};
    // I->addPolytope(6, H, h, OUTER);
    ;
};

int main() {
  /* to measure time */
  TicToc tt;

  /* construct grid for the state space */
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* grid node distance diameter */
  /* optimized values computed according to doi: 10.1109/CDC.2015.7403185 */
  double etaX[state_dim]={(25.0/362)*2*2*2,(3*M_PI/180/66)*2*2*2,(20.0/334)*2*2*2};
  /* lower bounds of the hyper rectangle */
  double lbX[state_dim]={58,-3*M_PI/180,0};
  /* upper bounds of the hyper rectangle */
  double ubX[state_dim]={83,0,20};
//  scots::UniformGrid ss(state_dim,s_lb,s_ub,s_eta);
//  std::cout << "Uniform grid details:" << std::endl;
//  ss.print_info();

  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
    double lbU[input_dim]={0,0};
  /* upper bounds of the hyper rectangle */
  double ubU[input_dim]={64000,8*M_PI/180};
  /* grid node distance diameter */
  double etaU[input_dim]={32000,(9.0/8.0*M_PI/180)};
    
    double tau = 0.25*1.5*1.5*1.5;
    int nSubInt = 5;
    
    double etaRatio[state_dim] = {2, 2, 2};
    double tauRatio = 1.5;
    int numAbs = 4;
    int p = 2;
    
    state_type x;
    input_type u;
    
    System aircraft(state_dim, lbX, ubX, etaX, tau,
                    input_dim, lbU, ubU, etaU,
                    etaRatio, tauRatio, nSubInt, numAbs);
    
    AdaptAbsReach abs("aircraft4A.log");
    abs.initialize(&aircraft, aircraftAddO, aircraftAddG);
    
    TicToc tt_tot;
    tt_tot.tic();
    abs.onTheFlyReach(p, aircraft_post, radius_post, x, u);
    clog << "------------------------------------Total time:";
    tt_tot.toc();
    
//  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
//  is.print_info();

  /* transition function of symbolic model */
//  scots::TransitionFunction tf;

  /* setup object to compute the transition function */
//  scots::Abstraction<state_type,input_type> abs(ss,is);
  /* measurement disturbances  */
//  state_type z={{0.0125,0.0025/180*M_PI,0.05}};
//  abs.set_measurement_error_bound(z);

//  std::cout << "Computing the transition function: " << std::endl;
//  tt.tic();
//  abs.compute_gb(tf,aircraft_post,radius_post);
//  tt.toc();
//  if(!getrusage(RUSAGE_SELF, &usage))
//    std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << std::endl;
//  std::cout << "Number of transitions: " << tf.get_no_transitions() << std::endl;

//  /* define target set */
//  auto target = [&s_eta, &z, &ss](const scots::abs_type& abs_state) {
//    state_type t_lb = {{63,-3*M_PI/180,0}};
//    state_type t_ub = {{75,0,2.5}};
//    state_type c_lb;
//    state_type c_ub;
//    /* center of cell associated with abs_state is stored in x */
//    state_type x;
//    ss.itox(abs_state,x);
//    /* hyper-interval of the quantizer symbol with perturbation */
//    for(int i=0; i<state_dim; i++) {
//      c_lb[i] = x[i]-s_eta[i]/2.0-z[i];
//      c_ub[i] = x[i]+s_eta[i]/2.0+z[i];
//    }
//    if( t_lb[0]<=c_lb[0] && c_ub[0]<=t_ub[0] &&
//        t_lb[1]<=c_lb[1] && c_ub[1]<=t_ub[1] &&
//        t_lb[2]<=c_lb[2] && c_ub[2]<=t_ub[2]) {
//      if(-0.91<=(x[0]*std::sin(x[1])-s_eta[0]/2.0-z[0]-(c_ub[0])*(s_eta[1]/2.0-z[1]))) {
//        return true;
//      }
//    }
//    return false;
//  };
//  /* write grid point IDs with uniform grid information to file */
//  write_to_file(ss,target,"target");
 
//  std::cout << "\nSynthesis: " << std::endl;
//  tt.tic();
//  scots::WinningDomain win=scots::solve_reachability_game(tf,target);
//  tt.toc();
//  std::cout << "Winning domain size: " << win.get_size() << std::endl;
//
//  std::cout << "\nWrite controller to controller.scs \n";
//  if(write_to_file(scots::StaticController(ss,is,std::move(win)),"controller"))
//    std::cout << "Done. \n";
//
//  return 1;
}
