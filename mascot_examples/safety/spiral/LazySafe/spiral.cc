#include <array>
#include <iostream>
#include <math.h>
#include <cmath>

#include "AdaptAbsSafe.hh"

#define _USE_MATH_DEFINES


using namespace std;
using namespace scots;
//using namespace helper;

/* state space dim */
#define dimX 2
#define dimU 1

/* disturbance */
const double w[2] = {0, 0};

/* angular speed of the spiral (fixed) */
const double omega = 1;

/* rate of change of the radius */
const double a_ = 0.1;

/* maximum radiun */
const double maxRadius = 2;

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* system ODE (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    /* cartesian to polar conversion */
    X_type y;
    y[0] = pow(pow(x[0], 2) + pow(x[1], 2), 0.5); // radius
    y[1] = atan2(x[1],x[0]); // angle
    auto sysODE = [](X_type &dydt, const X_type &y, const U_type &u) -> void {
        dydt[0] = -a_*y[0] + u[0]; /* radius in polar coordinate */
        dydt[1] = omega; /* angular velocity */
    };
    solver(sysODE, y, u);
//    if (x[0] > maxRadius) {
//        x[0] = maxRadius;
//    }
    if (y[0] < 0) {
        y[0] = -y[0];
        y[1] = y[1] + M_PI;
    }
    y[1] = fmod(y[1],(2*M_PI));
    /* polar to cartesian conversion */
    x[0] = y[0]*cos(y[1]);
    x[1] = y[0]*sin(y[1]);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    /* cartesian to polar conversion */
    X_type s;
    s[0] = pow(pow(r[0], 2) + pow(r[1], 2), 0.5);
    s[1] = atan2(r[1],r[0]);
    auto radODE = [](X_type &dsdt, const X_type &s, const U_type &u) -> void {
        double L[2][2];
        L[0][0] = -a_;
        L[0][1] = 0;
        L[1][0] = 0;
        L[1][1] = 0;

        dsdt[0] = L[0][0]*s[0] + L[0][1]*s[1] + w[0];
        dsdt[1] = L[1][0]*s[0] + L[1][1]*s[1] + w[1];
    };
    solver(radODE, s, u);
    /* polar to cartesian conversion */
    r[0] = s[0]*cos(s[1]);
    r[1] = s[0]*sin(s[1]);
};

//auto spiralAddS = [](SymbolicSet* S) -> void {
//    /* construct SymbolicSet for the safe set */
//    double H[4*dimX]={-1, 0,
//        1, 0,
//        0,-1,
//        0, 1};
//    /* add outer approximation of P={ x | H x<= h } form state space */
//    double h1[4] = {maxRadius,maxRadius,maxRadius,maxRadius};
//    /* initialize the safe set with the ss
//     * in order to obtain all the necessary information */
//    S->addPolytope(4, H, h1, INNER);
//    double h2[4] = {maxRadius,-0.775*maxRadius,0.1*maxRadius,0.1*maxRadius};
//    S->remPolytope(4, H, h2, OUTER);
//    double h3[4] = {0.725*maxRadius,-0.525*maxRadius,0.1*maxRadius,0.1*maxRadius};
//    S->remPolytope(4, H, h3, OUTER);
//    double h4[4] = {0.475*maxRadius,-0.275*maxRadius,0.1*maxRadius,0.1*maxRadius};
//    S->remPolytope(4, H, h4, OUTER);
//    double h5[4] = {0.225*maxRadius,-0.025*maxRadius,0.1*maxRadius,0.1*maxRadius};
//    S->remPolytope(4, H, h5, OUTER);
//};

auto spiralAddS = [](SymbolicSet* S) -> void {
    /* construct SymbolicSet for the safe set */
    double H[4*dimX]={-1, 0,
        1, 0,
        0,-1,
        0, 1};
    /* add outer approximation of P={ x | H x<= h } form state space */
    double h1[4] = {maxRadius,maxRadius,maxRadius,maxRadius};
    /* initialize the safe set with the ss
     * in order to obtain all the necessary information */
    S->addPolytope(4, H, h1, INNER);
    double h2[4] = {maxRadius,-0.8*maxRadius,0.025*maxRadius,0.025*maxRadius};
    S->remPolytope(4, H, h2, OUTER);
    double h3[4] = {0.7*maxRadius,-0.55*maxRadius,0.025*maxRadius,0.025*maxRadius};
    S->remPolytope(4, H, h3, OUTER);
    double h4[4] = {0.45*maxRadius,-0.3*maxRadius,0.025*maxRadius,0.025*maxRadius};
    S->remPolytope(4, H, h4, OUTER);
    double h5[4] = {0.2*maxRadius,-0.1*maxRadius,0.025*maxRadius,0.025*maxRadius};
    S->remPolytope(4, H, h5, OUTER);
};

int main() {

    double lbX[dimX]={-maxRadius, -maxRadius};
    double ubX[dimX]={maxRadius, maxRadius};
    double etaX[dimX]= {0.025*2*2, 0.025*2*2};

    double lbU[dimU]={-0.5};
    double ubU[dimU]= {0.5};
    double etaU[dimU]= {0.05};
    
    /* sampling time */
    const double tau = 0.05;
    /* number of intermediate steps in the ode solver */
    const int nint=5;
    
    double etaRatio[dimX] = {2, 2};
    double tauRatio = 1;
    
    int numAbs = 3;
    int readAbs = 0;
    
    X_type x;
    U_type u;
    
    System spiral(dimX, lbX, ubX, etaX, tau,
                  dimU, lbU, ubU, etaU,
                  etaRatio, tauRatio, nint, numAbs);
    
    AdaptAbsSafe abs("spiral_3A.log");
    abs.initialize(&spiral, spiralAddS);
    
    TicToc timer;
    timer.tic();
    abs.onTheFlySafe(sysNext, radNext, x, u);
    clog << "-----------------------------------------------Total time: " << timer.toc() << " seconds.\n";
    
    return 1;
}
