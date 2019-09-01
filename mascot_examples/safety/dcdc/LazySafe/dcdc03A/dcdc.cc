#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "AdaptAbsSafe.hh"

using namespace scots;

/* dimensions */
#define dimX 2
#define dimU 1

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

const double w[dimX] = {0.001, 0.001};

/* we integrate the simple ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, OdeSolver solver) -> void {
    auto sysODE = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
        const double r0 = 1.0 ;
        const double vs = 1.0 ;
        const double rl = 0.05 ;
        const double rc = rl / 10 ;
        const double xl = 3.0 ;
        const double xc = 70.0 ;

        const double b[2]={vs/xl, 0};

        double a[2][2];
        if(u[0]==0.5) {
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
auto radNext = [](X_type &r, U_type &u, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        const double r0 = 1.0 ;
        const double rl = 0.05 ;
        const double rc = rl / 10 ;
        const double xl = 3.0 ;
        const double xc = 70.0 ;

        double a[2][2];
        if(u[0]==0.5) {
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

auto dcdcAddO = [](SymbolicSet* O) -> void {
    ;
};

auto dcdcAddS = [](SymbolicSet* S) -> void {
    double H[4*2] = {-1,  0,
                         1,  0,
                         0, -1,
                         0,  1};
    double h[4] = {-1.15, 1.55, -5.45, 5.85};
    S->addPolytope(4, H, h, INNER);
};

int main() {

    double lbX[dimX]  = {1.15, 5.45};
    double ubX[dimX]  = {1.55, 5.85};

    /* the system dynamics has two modes which correspond to two distinct abstract control inputs (in our case, they are 0.5, 1.5) */
    double lbU[dimU]  = {0};
    double ubU[dimU]  = {2};
    double etaU[dimU] = {1};

    int nint = 5;

//    double etaX[dimX]= {(pow(2,2)*2/4e3), (pow(2,2)*2/4e3)};
    double etaX[dimX]= {(pow(2,2)*2*2*2/4e3), (pow(2,2)*2*2*2/4e3)};
    //double tau = pow(2, 2)*0.0625;
    double tau = 0.5;
    int numAbs = 3;

    double etaRatio[dimX] = {2, 2};
    double tauRatio = 1;
    int verbose = 1;

    X_type x;
    U_type u;


    System dcdc(dimX, lbX, ubX, etaX, tau,
                dimU, lbU, ubU, etaU,
                etaRatio, tauRatio, nint, numAbs);
    AdaptAbsSafe abs("dcdc3A.log", verbose);
    abs.initialize(&dcdc, dcdcAddS);

    TicToc timer;
    timer.tic();
    abs.onTheFlySafeNoAux(sysNext, radNext, x, u);
    clog << "-----------------------------------------------Total time: " << timer.toc() << " seconds.\n";
}
