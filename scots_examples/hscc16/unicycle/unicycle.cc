/*
 * unicycle.cc
 *
 *  created on: 21.01.2016
 *      author: rungger
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <array>
#include <iostream>

#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"


/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/* sampling time */
const double tau = 0.3;
/* number of intermediate steps in the ode solver */
const int nint=5;
OdeSolver ode_solver(sDIM,nint,tau);

/* we integrate the unicycle ode by 0.3 sec (the result is stored in x)  */
auto  unicycle_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the unicycle */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) -> void {
      xx[0] = u[0]*std::cos(x[2]);
      xx[1] = u[0]*std::sin(x[2]);
      xx[2] = u[1];
  };
  ode_solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) -> void {
    r[0] = r[0]+r[2]*std::abs(u[0])*0.3;
    r[1] = r[1]+r[2]*std::abs(u[0])*0.3;
};


/* forward declaration of the functions to setup the state space 
 * and input space of the unicycle example */
scots::SymbolicSet unicycleCreateStateSpace(Cudd &mgr);
scots::SymbolicSet unicycleCreateInputSpace(Cudd &mgr);


int main() {
  /* to measure time */
  TicToc tt;
  /* there is one unique manager to organize the bdd variables */
  Cudd mgr;

  /****************************************************************************/
  /* construct SymbolicSet for the state space */
  /****************************************************************************/
  scots::SymbolicSet ss=unicycleCreateStateSpace(mgr);
  ss.writeToFile("unicycle_ss.bdd");
  /* write SymbolicSet of obstacles to unicycle_obst.bdd */
  ss.complement();
  ss.writeToFile("unicycle_obst.bdd");
  ss.complement();
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* the target set */
  /****************************************************************************/
  /* first make a copy of the state space so that we obtain the grid
   * information in the new symbolic set */
  scots::SymbolicSet ts = ss;
  /* define the target set as a symbolic set */
  double H[9]={ 2, 0, 0,
                0, 1, 0,
                0, 0, .1};
  /* compute inner approximation of P={ x | H x<= h1 }  */
  double c[3] = {9.5,0.6,0};
  ts.addEllipsoid(H,c, scots::INNER);
  ts.writeToFile("unicycle_target.bdd");
    std::cout << "ts:\n";
    ts.printInfo(1);



  /****************************************************************************/
  /* construct SymbolicSet for the input space */
  /****************************************************************************/
  scots::SymbolicSet is=unicycleCreateInputSpace(mgr);
  std::cout << std::endl << "Input space details:" << std::endl;
  is.printInfo(1);

  /****************************************************************************/
  /* setup class for symbolic model computation */
  /****************************************************************************/
  /* first create SymbolicSet of post variables 
   * by copying the SymbolicSet of the state space and assigning new BDD IDs */
  scots::SymbolicSet sspost(ss,1);
  /* instantiate the SymbolicModel */
  scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
  /* compute the transition relation */
  tt.tic();
  abstraction.computeTransitionRelation(unicycle_post, radius_post);
  std::cout << std::endl;
  tt.toc();
  /* get the number of elements in the transition relation */
  std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;


  /****************************************************************************/
  /* we continue with the controller synthesis */
  /****************************************************************************/
  int verbose=1;
  /* we setup a fixed point object to compute reachabilty controller */
  scots::FixedPoint fp(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD T = ts.getSymbolicSet();
  tt.tic();
  BDD C = fp.reach(T,verbose);
  tt.toc();

  /****************************************************************************/
  /* last we store the controller as a SymbolicSet 
   * the underlying uniform grid is given by the Cartesian product of 
   * the uniform gird of the space and uniform gird of the input space */
  /****************************************************************************/
  scots::SymbolicSet controller(ss,is);
  controller.setSymbolicSet(C);
  controller.writeToFile("unicycle_controller.bdd");

  scots::SymbolicSet tr = abstraction.getTransitionRelation();
  tr.writeToFile("unicycle_transition_relation.bdd");

    std::cout << "tr:\n";
    tr.printInfo(1);
    std::cout << "controller:\n";
    controller.printInfo(1);
  
  return 1;
}

scots::SymbolicSet unicycleCreateStateSpace(Cudd &mgr) {

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={0,0,-M_PI-0.4};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={10,10,M_PI+0.4}; 
  /* grid node distance diameter */
  double eta[sDIM]={.2,.2,.1};   



  /* eta is added to the bound so as to ensure that the whole
   * [0,10]x[0,10]x[-pi-eta,pi+eta] is covered by the cells */

  scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);

  /* add the grid points to the SymbolicSet ss */
  ss.addGridPoints();
  /* remove the obstacles from the state space */
  /* the obstacles are defined as polytopes */
  /* define H* x <= h */
  double H[4*sDIM]={-1, 0, 0,
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
  /* remove outer approximation of P={ x | H x<= h1 } form state space */
  double h1[4] = {-1,1.2,-0, 9};
  ss.remPolytope(4,H,h1, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h2 } form state space */
  double h2[4] = {-2.2,2.4,-0,5};
  ss.remPolytope(4,H,h2, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h3 } form state space */
  double h3[4] = {-2.2,2.4,-6,10};
  ss.remPolytope(4,H,h3, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h4 } form state space */
  double h4[4] = {-3.4,3.6,-0,9};
  ss.remPolytope(4,H,h4, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h5 } form state space */
  double h5[4] = {-4.6 ,4.8,-1,10};
  ss.remPolytope(4,H,h5, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h6 } form state space */
  double h6[4] = {-5.8,6,-0,6};
  ss.remPolytope(4,H,h6, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h7 } form state space */
  double h7[4] = {-5.8,6,-7,10};
  ss.remPolytope(4,H,h7, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h8 } form state space */
  double h8[4] = {-7,7.2,-1,10};
  ss.remPolytope(4,H,h8, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h9 } form state space */
  double h9[4] = {-8.2,8.4,-0,8.5};
  ss.remPolytope(4,H,h9, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h10 } form state space */
  double h10[4] = {-8.4,9.3,-8.3,8.5};
  ss.remPolytope(4,H,h10, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h11 } form state space */
  double h11[4] = {-9.3,10,-7.1,7.3};
  ss.remPolytope(4,H,h11, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h12 } form state space */
  double h12[4] = {-8.4,9.3,-5.9,6.1};
  ss.remPolytope(4,H,h12, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h13 } form state space */
  double h13[4] = {-9.3,10 ,-4.7,4.9};
  ss.remPolytope(4,H,h13, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h14 } form state space */
  double h14[4] = {-8.4,9.3,-3.5,3.7};
  ss.remPolytope(4,H,h14, scots::OUTER);
  /* remove outer approximation of P={ x | H x<= h15 } form state space */
  double h15[4] = {-9.3,10 ,-2.3,2.5};
  ss.remPolytope(4,H,h15, scots::OUTER);


 return ss;
}

scots::SymbolicSet unicycleCreateInputSpace(Cudd &mgr) {

  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={-1,-1.5};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={1,1.5}; 
  /* grid node distance diameter */
  double eta[sDIM]={.3,.2};   

  scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
  is.addGridPoints();

  return is;
}

