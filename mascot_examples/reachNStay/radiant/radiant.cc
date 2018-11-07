#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsReach.hh"
#include "AdaptAbsSafe.hh"

using namespace scots;

/* dimensions */
#define dimX 3
#define dimU 1

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

const double w[dimX] = {0, 0}; // disturbance set to 0 for comparison with cosyma

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, OdeSolver solver) -> void {
    auto sysODE = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
//        const double r0 = 1.0 ;
//        const double vs = 1.0 ;
//        const double rl = 0.05 ;
//        const double rc = rl / 10 ;
//        const double xl = 3.0 ;
//        const double xc = 70.0 ;
//
//        const double b[2]={vs/xl, 0};
//        -0.0413    0.0106    0.0080
//        0.4377   -0.4869    0.0260
//        0.4377    0.0346   -0.4955
//
//        -0.0186    0.0106    0.0080
//        0.4377   -0.4869    0.0260
//        0.4377    0.0346   -0.4955
        
        double a[3][3];
        double b[3];
        if(u[0]==0.5) {
            a[0][0] = -0.0413*0.001;
            a[0][1] = 0.0106*0.001;
            a[0][2] = 0.0080*0.001;
            a[1][0] = 0.4377*0.001;
            a[1][1] = -0.4869*0.001;
            a[1][2] = 0.0260*0.001;
            a[2][0] = 0.4377*0.001;
            a[2][1] = 0.0346*0.001;
            a[2][2] = -0.4955*0.001;
            
            b[0] = 0.0004;
            b[1] = 0.0010;
            b[2] = 0.0011;
        } else {
            a[0][0] = -0.0186*0.001;
            a[0][1] = 0.0106*0.001;
            a[0][2] = 0.0080*0.001;
            a[1][0] = 0.4377*0.001;
            a[1][1] = -0.4869*0.001;
            a[1][2] = 0.0260*0.001;
            a[2][0] = 0.4377*0.001;
            a[2][1] = 0.0346*0.001;
            a[2][2] = -0.4955*0.001;
            
            b[0] = 0;
            b[1] = 0.0010;
            b[2] = 0.0011;
        }
        dxdt[0] = a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2] + b[0];
        dxdt[1] = a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2] + b[1];
        dxdt[2] = a[2][0]*x[0] + a[2][1]*x[1] + a[2][2]*x[2] + b[2];
    };
    solver(sysODE, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        double a[3][3];
        double b[3];
        if(u[0]==0.5) {
            a[0][0] = -0.0413*0.001;
            a[0][1] = 0.0106*0.001;
            a[0][2] = 0.0080*0.001;
            a[1][0] = 0.4377*0.001;
            a[1][1] = -0.4869*0.001;
            a[1][2] = 0.0260*0.001;
            a[2][0] = 0.4377*0.001;
            a[2][1] = 0.0346*0.001;
            a[2][2] = -0.4955*0.001;
        } else {
            a[0][0] = -0.0186*0.001;
            a[0][1] = 0.0106*0.001;
            a[0][2] = 0.0080*0.001;
            a[1][0] = 0.4377*0.001;
            a[1][1] = -0.4869*0.001;
            a[1][2] = 0.0260*0.001;
            a[2][0] = 0.4377*0.001;
            a[2][1] = 0.0346*0.001;
            a[2][2] = -0.4955*0.001;
        }
        drdt[0] = a[0][0]*r[0] + a[0][1]*r[1] + a[0][2]*r[2];
        drdt[1] = a[1][0]*r[0] + a[1][1]*r[1] + a[1][2]*r[2];
        drdt[2] = a[2][0]*r[0] + a[2][1]*r[1] + a[2][2]*r[2];
    };
    solver(radODE, r, u);
};

auto dcdcAddO = [](SymbolicSet* O) -> void {
    ;
};

auto dcdcAddS = [](SymbolicSet* S) -> void {
    double H[6*3] = {-1,  0, 0,
                         1,  0, 0,
                         0, -1, 0,
                         0,  1, 0,
                         0,  0, -1,
                         0,  0, 1};
    double h[6] = {-21, 27, -22, 25, -22, 25};
    S->addPolytope(6, H, h, INNER);
};

auto dcdcAddG = [](SymbolicSet* G) -> void {
    double H[4*2] = {-1,  0,
        1,  0,
        0, -1,
        0,  1};
    double h[4] = {-1.1, 1.6, -5.4, 5.9};
    G->addPolytope(4, H, h, INNER);
};

int main() {
    
    double lbX[dimX]  = {20, 20, 20};
    double ubX[dimX]  = {28, 28, 28};
    
    /* the system dynamics has two modes which correspond to two distinct abstract control inputs (in our case, they are 0.5, 1.5) */
    double lbU[dimU]  = {0};
    double ubU[dimU]  = {2};
    double etaU[dimU] = {1};
    
    int nint = 5;
    
    double etaX[dimX]= {pow(10,-3), pow(10,-3), pow(10,-3)};
    //double tau = pow(2, 2)*0.0625;
    double tau = 2;
    int numAbs = 3;
    
    double etaRatio[dimX] = {2, 2};
    double tauRatio = 2;
    
    X_type x;
    U_type u;
    
    
    System dcdc(dimX, lbX, ubX, etaX, tau,
                dimU, lbU, ubU, etaU,
                etaRatio, tauRatio, nint, numAbs);
    int verbose = 0;
    int p = 3;
    AdaptAbsReach abs("dcdc3A.log", verbose);
    abs.initialize(&dcdc, dcdcAddO, dcdcAddG);
    
    TicToc timer;
    timer.tic();
    abs.onTheFlyReach(p,sysNext, radNext, x, u);
    clog << "-----------------------------------------------Total time: " << timer.toc() << " seconds.\n";
}

