#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "Adaptive.hh"
#include "Compare.hh"

using namespace scots;

/* dimensions */
#define dimX 2
#define dimU 1

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

const double w[dimX] = {0.001, 0.001};

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto sysODE = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
        const double r0 = 1.0 ;
        const double vs = 1.0 ;
        const double rl = 0.05 ;
        const double rc = rl / 10 ;
        const double xl = 3.0 ;
        const double xc = 70.0 ;

        const double b[2]={vs/xl, 0};

        double a[2][2];
        if(u[0]==1) {
            a[0][0] = -rl / xl;
            a[0][1] = 0;
            a[1][0] = 0;
            a[1][1] = (-1 / xc) * (1 / (r0 + rc));
        } else {
            a[0][0] = (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc))) ;
            a[0][1] =  ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
            a[1][0] = 5 * (r0 / (r0 + rc)) * (1 / xc);
            a[1][1] =(-1 / xc) * (1 / (r0 + rc)) ;
        }
        dxdt[0] = a[0][0]*x[0]+a[0][1]*x[1] + b[0];
        dxdt[1] = a[1][0]*x[0]+a[1][1]*x[1] + b[1];
    };
    solver(sysODE, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        const double r0 = 1.0 ;
        const double rl = 0.05 ;
        const double rc = rl / 10 ;
        const double xl = 3.0 ;
        const double xc = 70.0 ;

        double a[2][2];
        if(u[0]==1) {
            a[0][0] = -rl / xl;
            a[0][1] = 0;
            a[1][0] = 0;
            a[1][1] = (-1 / xc) * (1 / (r0 + rc));
        } else {
            a[0][0] = (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc))) ;
            a[0][1] =  ((1 / xl) * (r0 / (r0 + rc))) / 5 ;
            a[1][0] = 5 * (r0 / (r0 + rc)) * (1 / xc);
            a[1][1] =(-1 / xc) * (1 / (r0 + rc)) ;
        }
        drdt[0] = a[0][0]*r[0]+a[0][1]*r[1] + w[0];
        drdt[1] = a[1][0]*r[0]+a[1][1]*r[1] + w[1];
    };
    solver(radODE, r, u);
};

auto dcdcAddS = [](SymbolicSet* S) -> void {
    double H[4*2] = {-1,  0,
                         1,  0,
                         0, -1,
                         0,  1};
    double h[4] = {-1.15, 1.55, -5.45, 5.85};
    S->addPolytope(4, H, h, INNER);
};

void sub(double* lbX, double* ubX, double* etaX, double tau, double* lbU, double* ubU, double* etaU,
         double* etaRatio, double tauRatio, int numAbs, int nint, int readAb) {
    Compare<X_type, U_type> comp(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readAb, "scots.txt");
    comp.initializeSafe(dcdcAddS);
    comp.computeAbstractions(sysNext, radNext);
    comp.safeSCOTS();
}

int main() {

    double lbX[dimX]  = {1.15, 5.45};
    double ubX[dimX]  = {1.55, 5.85};

    double lbU[dimU]  = {1};
    double ubU[dimU]  = {2};
    double etaU[dimU] = {1};

    int nint = 5;

    double etaX[dimX]= {2/4e3*3*3, 2/4e3*3*3};
    double tau = 0.5;

    double etaRatio[dimX] = {3, 3};
    double tauRatio = 1;

    int numAbs = 2;
    int readXX = 0; // if X, U have changed, needs to be 0.
    int readAbs = 0; // if above or dynamics have changed, needs to be 0.

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readXX, readAbs, "adaptive.txt");
    abs.initializeSafe(dcdcAddS);
    abs.computeAbstractions(sysNext, radNext);
    abs.safe();

    sub(lbX, ubX, etaX, tau, lbU, ubU, etaU, etaRatio, tauRatio, numAbs, nint, readAbs);
}
