#include <array>
#include <iostream>
#include <cmath>
#define _USE_MATH_DEFINES

#include "Adaptive.hh"


using namespace std;
using namespace scots;

/* state space dim */
int dimX = 4;
int dimU = 1;

/* data types for the ode solver */
typedef std::array<double,4> X_type;
typedef std::array<double,1> U_type;

/* mechanical parameters */
const double l1 = 0.323; // length
const double l2 = 0.480; // [m]
const double a1 = 0.2145; // distance to center of gravity
const double a2 = 0.223; // [m]
const double m1 = 0.853; // mass
const double m2 = 0.510; // [kg]
const double J1 = 0.0126; // moment of inertia
const double J2 = 0.0185; // [N m s2]
const double d1 = 0.005; // friction constant
const double d2 = 0.005; // [N m s]
const double g = 9.81; // standard gravity [m s^-2]

/* we integrate the pendulum ode by 0.3 sec (the result is stored in x)  */
auto sysNext = [](X_type &x, U_type &u, double tau, OdeSolver solver) -> void {
    auto sysODE = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
        // u[0] = u[0]
        // x[0] = x[0]
        // x[1] = x[1]
        dxdt[0] = x[2]; // x[2]
        dxdt[1] = x[3]; // x[3]
        dxdt[2] = (J2*d2*x[3] - J2*d2*x[2] - J2*d1*x[2] - pow(a2,2)*d1*x[2]*m2 - pow(a2,2)*d2*x[2]*m2 + pow(a2,2)*d2*x[3]*m2 - pow(a2,3)*pow(x[3],2)*l1*pow(m2,2)*sin(x[0] - x[1]) + (pow(a2,2)*u[0]*l1*pow(m2,2)*cos(x[0]))/2 + (pow(a2,2)*g*l1*pow(m2,2)*sin(x[0]))/2 + J2*a1*u[0]*m1*cos(x[0]) + J2*u[0]*l1*m2*cos(x[0]) - (pow(a2,2)*u[0]*l1*pow(m2,2)*cos(x[0] - 2*x[1]))/2 + J2*a1*g*m1*sin(x[0]) + J2*g*l1*m2*sin(x[0]) + (pow(a2,2)*g*l1*pow(m2,2)*sin(x[0] - 2*x[1]))/2 - (pow(a2,2)*pow(x[2],2)*pow(l1,2)*pow(m2,2)*sin(2*x[0] - 2*x[1]))/2 - a2*d2*x[2]*l1*m2*cos(x[0] - x[1]) + a2*d2*x[3]*l1*m2*cos(x[0] - x[1]) + a1*pow(a2,2)*u[0]*m1*m2*cos(x[0]) + a1*pow(a2,2)*g*m1*m2*sin(x[0]) - J2*a2*pow(x[3],2)*l1*m2*sin(x[0] - x[1]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(x[0] - x[1]), 2) + pow(a1,2)*pow(a2,2)*m1*m2);
        dxdt[3] = (J1*d2*x[2] - J1*d2*x[3] + pow(a1,2)*d2*x[2]*m1 - pow(a1,2)*d2*x[3]*m1 + d2*x[2]*pow(l1,2)*m2 - d2*x[3]*pow(l1,2)*m2 + a2*pow(x[2],2)*pow(l1,3)*pow(m2,2)*sin(x[0] - x[1]) + a2*u[0]*pow(l1,2)*pow(m2,2)*cos(x[1]) + a2*g*pow(l1,2)*pow(m2,2)*sin(x[1]) + J1*a2*u[0]*m2*cos(x[1]) + J1*a2*g*m2*sin(x[1]) + (pow(a2,2)*pow(x[3],2)*pow(l1,2)*pow(m2,2)*sin(2*x[0] - 2*x[1]))/2 + a2*d1*x[2]*l1*m2*cos(x[0] - x[1]) + a2*d2*x[2]*l1*m2*cos(x[0] - x[1]) - a2*d2*x[3]*l1*m2*cos(x[0] - x[1]) + pow(a1,2)*a2*u[0]*m1*m2*cos(x[1]) - a2*u[0]*pow(l1,2)*pow(m2,2)*cos(x[0] - x[1])*cos(x[0]) + pow(a1,2)*a2*g*m1*m2*sin(x[1]) - a2*g*pow(l1,2)*pow(m2,2)*cos(x[0] - x[1])*sin(x[0]) + J1*a2*pow(x[2],2)*l1*m2*sin(x[0] - x[1]) + pow(a1,2)*a2*pow(x[2],2)*l1*m1*m2*sin(x[0] - x[1]) - a1*a2*u[0]*l1*m1*m2*cos(x[0] - x[1])*cos(x[0]) - a1*a2*g*l1*m1*m2*cos(x[0] - x[1])*sin(x[0]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(x[0] - x[1]), 2) + pow(a1,2)*pow(a2,2)*m1*m2);
    };
    solver(sysODE, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radNext = [](X_type &r, U_type &u, double tau, OdeSolver solver) -> void {
    auto radODE = [](X_type &drdt, const X_type &r, const U_type &u) -> void {
        double L[4][4];
        L[0][0] = 0;
        L[0][1] = 0;
        L[0][2] = 1;
        L[0][3] = 0;
        L[1][0] = 0;
        L[1][1] = 0;
        L[1][2] = 0;
        L[1][3] = 1;
        L[2][0] = abs( - (pow(a2,3)*pow(r[3],2)*l1*pow(m2,2)*cos(r[0] - r[1]) - (pow(a2,2)*g*l1*pow(m2,2)*cos(r[0]))/2 + (pow(a2,2)*u[0]*l1*pow(m2,2)*sin(r[0]))/2 - J2*a1*g*m1*cos(r[0]) - J2*g*l1*m2*cos(r[0]) + J2*a1*u[0]*m1*sin(r[0]) - (pow(a2,2)*g*l1*pow(m2,2)*cos(r[0] - 2*r[1]))/2 + J2*u[0]*l1*m2*sin(r[0]) - (pow(a2,2)*u[0]*l1*pow(m2,2)*sin(r[0] - 2*r[1]))/2 + pow(a2,2)*pow(r[2],2)*pow(l1,2)*pow(m2,2)*cos(2*r[0] - 2*r[1]) - a2*d2*r[2]*l1*m2*sin(r[0] - r[1]) + a2*d2*r[3]*l1*m2*sin(r[0] - r[1]) - a1*pow(a2,2)*g*m1*m2*cos(r[0]) + a1*pow(a2,2)*u[0]*m1*m2*sin(r[0]) + J2*a2*pow(r[3],2)*l1*m2*cos(r[0] - r[1]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2) - (2*pow(a2,2)*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*sin(r[0] - r[1])*(J2*d2*r[3] - J2*d2*r[2] - J2*d1*r[2] - pow(a2,2)*d1*r[2]*m2 - pow(a2,2)*d2*r[2]*m2 + pow(a2,2)*d2*r[3]*m2 - pow(a2,3)*pow(r[3],2)*l1*pow(m2,2)*sin(r[0] - r[1]) + (pow(a2,2)*u[0]*l1*pow(m2,2)*cos(r[0]))/2 + (pow(a2,2)*g*l1*pow(m2,2)*sin(r[0]))/2 + J2*a1*u[0]*m1*cos(r[0]) + J2*u[0]*l1*m2*cos(r[0]) - (pow(a2,2)*u[0]*l1*pow(m2,2)*cos(r[0] - 2*r[1]))/2 + J2*a1*g*m1*sin(r[0]) + J2*g*l1*m2*sin(r[0]) + (pow(a2,2)*g*l1*pow(m2,2)*sin(r[0] - 2*r[1]))/2 - (pow(a2,2)*pow(r[2],2)*pow(l1,2)*pow(m2,2)*sin(2*r[0] - 2*r[1]))/2 - a2*d2*r[2]*l1*m2*cos(r[0] - r[1]) + a2*d2*r[3]*l1*m2*cos(r[0] - r[1]) + a1*pow(a2,2)*u[0]*m1*m2*cos(r[0]) + a1*pow(a2,2)*g*m1*m2*sin(r[0]) - J2*a2*pow(r[3],2)*l1*m2*sin(r[0] - r[1])))/pow((m1*pow(a1,2)*pow(a2,2)*m2 + J2*m1*pow(a1,2) - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a2,2)*pow(l1,2)*pow(m2,2) + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 + J1*J2),2) );
        L[2][1] = abs( (pow(a2,3)*pow(r[3],2)*l1*pow(m2,2)*cos(r[0] - r[1]) - pow(a2,2)*g*l1*pow(m2,2)*cos(r[0] - 2*r[1]) - pow(a2,2)*u[0]*l1*pow(m2,2)*sin(r[0] - 2*r[1]) + pow(a2,2)*pow(r[2],2)*pow(l1,2)*pow(m2,2)*cos(2*r[0] - 2*r[1]) - a2*d2*r[2]*l1*m2*sin(r[0] - r[1]) + a2*d2*r[3]*l1*m2*sin(r[0] - r[1]) + J2*a2*pow(r[3],2)*l1*m2*cos(r[0] - r[1]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2) + (2*pow(a2,2)*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*sin(r[0] - r[1])*(J2*d2*r[3] - J2*d2*r[2] - J2*d1*r[2] - pow(a2,2)*d1*r[2]*m2 - pow(a2,2)*d2*r[2]*m2 + pow(a2,2)*d2*r[3]*m2 - pow(a2,3)*pow(r[3],2)*l1*pow(m2,2)*sin(r[0] - r[1]) + (pow(a2,2)*u[0]*l1*pow(m2,2)*cos(r[0]))/2 + (pow(a2,2)*g*l1*pow(m2,2)*sin(r[0]))/2 + J2*a1*u[0]*m1*cos(r[0]) + J2*u[0]*l1*m2*cos(r[0]) - (pow(a2,2)*u[0]*l1*pow(m2,2)*cos(r[0] - 2*r[1]))/2 + J2*a1*g*m1*sin(r[0]) + J2*g*l1*m2*sin(r[0]) + (pow(a2,2)*g*l1*pow(m2,2)*sin(r[0] - 2*r[1]))/2 - (pow(a2,2)*pow(r[2],2)*pow(l1,2)*pow(m2,2)*sin(2*r[0] - 2*r[1]))/2 - a2*d2*r[2]*l1*m2*cos(r[0] - r[1]) + a2*d2*r[3]*l1*m2*cos(r[0] - r[1]) + a1*pow(a2,2)*u[0]*m1*m2*cos(r[0]) + a1*pow(a2,2)*g*m1*m2*sin(r[0]) - J2*a2*pow(r[3],2)*l1*m2*sin(r[0] - r[1])))/pow((m1*pow(a1,2)*pow(a2,2)*m2 + J2*m1*pow(a1,2) - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a2,2)*pow(l1,2)*pow(m2,2) + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 + J1*J2),2) );
        L[2][2] = -(J2*d1 + J2*d2 + pow(a2,2)*d1*m2 + pow(a2,2)*d2*m2 + a2*d2*l1*m2*cos(r[0] - r[1]) + pow(a2,2)*r[2]*pow(l1,2)*pow(m2,2)*sin(2*r[0] - 2*r[1]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2);
        L[2][3] = abs( (J2*d2 + pow(a2,2)*d2*m2 + a2*d2*l1*m2*cos(r[0] - r[1]) - 2*pow(a2,3)*r[3]*l1*pow(m2,2)*sin(r[0] - r[1]) - 2*J2*a2*r[3]*l1*m2*sin(r[0] - r[1]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2) );
        L[3][0] = abs( (a2*pow(r[2],2)*pow(l1,3)*pow(m2,2)*cos(r[0] - r[1]) + pow(a2,2)*pow(r[3],2)*pow(l1,2)*pow(m2,2)*cos(2*r[0] - 2*r[1]) - a2*d1*r[2]*l1*m2*sin(r[0] - r[1]) - a2*d2*r[2]*l1*m2*sin(r[0] - r[1]) + a2*d2*r[3]*l1*m2*sin(r[0] - r[1]) - a2*g*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*cos(r[0]) + a2*u[0]*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*sin(r[0]) + a2*u[0]*pow(l1,2)*pow(m2,2)*sin(r[0] - r[1])*cos(r[0]) + J1*a2*pow(r[2],2)*l1*m2*cos(r[0] - r[1]) + a2*g*pow(l1,2)*pow(m2,2)*sin(r[0] - r[1])*sin(r[0]) + pow(a1,2)*a2*pow(r[2],2)*l1*m1*m2*cos(r[0] - r[1]) - a1*a2*g*l1*m1*m2*cos(r[0] - r[1])*cos(r[0]) + a1*a2*u[0]*l1*m1*m2*cos(r[0] - r[1])*sin(r[0]) + a1*a2*u[0]*l1*m1*m2*sin(r[0] - r[1])*cos(r[0]) + a1*a2*g*l1*m1*m2*sin(r[0] - r[1])*sin(r[0]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2) - (2*pow(a2,2)*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*sin(r[0] - r[1])*(J1*d2*r[2] - J1*d2*r[3] + pow(a1,2)*d2*r[2]*m1 - pow(a1,2)*d2*r[3]*m1 + d2*r[2]*pow(l1,2)*m2 - d2*r[3]*pow(l1,2)*m2 + a2*pow(r[2],2)*pow(l1,3)*pow(m2,2)*sin(r[0] - r[1]) + a2*u[0]*pow(l1,2)*pow(m2,2)*cos(r[1]) + a2*g*pow(l1,2)*pow(m2,2)*sin(r[1]) + J1*a2*u[0]*m2*cos(r[1]) + J1*a2*g*m2*sin(r[1]) + (pow(a2,2)*pow(r[3],2)*pow(l1,2)*pow(m2,2)*sin(2*r[0] - 2*r[1]))/2 + a2*d1*r[2]*l1*m2*cos(r[0] - r[1]) + a2*d2*r[2]*l1*m2*cos(r[0] - r[1]) - a2*d2*r[3]*l1*m2*cos(r[0] - r[1]) + pow(a1,2)*a2*u[0]*m1*m2*cos(r[1]) - a2*u[0]*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*cos(r[0]) + pow(a1,2)*a2*g*m1*m2*sin(r[1]) - a2*g*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*sin(r[0]) + J1*a2*pow(r[2],2)*l1*m2*sin(r[0] - r[1]) + pow(a1,2)*a2*pow(r[2],2)*l1*m1*m2*sin(r[0] - r[1]) - a1*a2*u[0]*l1*m1*m2*cos(r[0] - r[1])*cos(r[0]) - a1*a2*g*l1*m1*m2*cos(r[0] - r[1])*sin(r[0])))/pow((m1*pow(a1,2)*pow(a2,2)*m2 + J2*m1*pow(a1,2) - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a2,2)*pow(l1,2)*pow(m2,2) + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 + J1*J2),2) );
        L[3][1] = abs( (2*pow(a2,2)*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*sin(r[0] - r[1])*(J1*d2*r[2] - J1*d2*r[3] + pow(a1,2)*d2*r[2]*m1 - pow(a1,2)*d2*r[3]*m1 + d2*r[2]*pow(l1,2)*m2 - d2*r[3]*pow(l1,2)*m2 + a2*pow(r[2],2)*pow(l1,3)*pow(m2,2)*sin(r[0] - r[1]) + a2*u[0]*pow(l1,2)*pow(m2,2)*cos(r[1]) + a2*g*pow(l1,2)*pow(m2,2)*sin(r[1]) + J1*a2*u[0]*m2*cos(r[1]) + J1*a2*g*m2*sin(r[1]) + (pow(a2,2)*pow(r[3],2)*pow(l1,2)*pow(m2,2)*sin(2*r[0] - 2*r[1]))/2 + a2*d1*r[2]*l1*m2*cos(r[0] - r[1]) + a2*d2*r[2]*l1*m2*cos(r[0] - r[1]) - a2*d2*r[3]*l1*m2*cos(r[0] - r[1]) + pow(a1,2)*a2*u[0]*m1*m2*cos(r[1]) - a2*u[0]*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*cos(r[0]) + pow(a1,2)*a2*g*m1*m2*sin(r[1]) - a2*g*pow(l1,2)*pow(m2,2)*cos(r[0] - r[1])*sin(r[0]) + J1*a2*pow(r[2],2)*l1*m2*sin(r[0] - r[1]) + pow(a1,2)*a2*pow(r[2],2)*l1*m1*m2*sin(r[0] - r[1]) - a1*a2*u[0]*l1*m1*m2*cos(r[0] - r[1])*cos(r[0]) - a1*a2*g*l1*m1*m2*cos(r[0] - r[1])*sin(r[0])))/pow((m1*pow(a1,2)*pow(a2,2)*m2 + J2*m1*pow(a1,2) - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a2,2)*pow(l1,2)*pow(m2,2) + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 + J1*J2),2) - (a2*pow(r[2],2)*pow(l1,3)*pow(m2,2)*cos(r[0] - r[1]) - a2*g*pow(l1,2)*pow(m2,2)*cos(r[1]) + a2*u[0]*pow(l1,2)*pow(m2,2)*sin(r[1]) - J1*a2*g*m2*cos(r[1]) + J1*a2*u[0]*m2*sin(r[1]) + pow(a2,2)*pow(r[3],2)*pow(l1,2)*pow(m2,2)*cos(2*r[0] - 2*r[1]) - a2*d1*r[2]*l1*m2*sin(r[0] - r[1]) - a2*d2*r[2]*l1*m2*sin(r[0] - r[1]) + a2*d2*r[3]*l1*m2*sin(r[0] - r[1]) - pow(a1,2)*a2*g*m1*m2*cos(r[1]) + pow(a1,2)*a2*u[0]*m1*m2*sin(r[1]) + a2*u[0]*pow(l1,2)*pow(m2,2)*sin(r[0] - r[1])*cos(r[0]) + J1*a2*pow(r[2],2)*l1*m2*cos(r[0] - r[1]) + a2*g*pow(l1,2)*pow(m2,2)*sin(r[0] - r[1])*sin(r[0]) + pow(a1,2)*a2*pow(r[2],2)*l1*m1*m2*cos(r[0] - r[1]) + a1*a2*u[0]*l1*m1*m2*sin(r[0] - r[1])*cos(r[0]) + a1*a2*g*l1*m1*m2*sin(r[0] - r[1])*sin(r[0]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2) );
        L[3][2] = abs( (J1*d2 + pow(a1,2)*d2*m1 + d2*pow(l1,2)*m2 + a2*d1*l1*m2*cos(r[0] - r[1]) + a2*d2*l1*m2*cos(r[0] - r[1]) + 2*a2*r[2]*pow(l1,3)*pow(m2,2)*sin(r[0] - r[1]) + 2*J1*a2*r[2]*l1*m2*sin(r[0] - r[1]) + 2*pow(a1,2)*a2*r[2]*l1*m1*m2*sin(r[0] - r[1]))/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2) );
        L[3][3] = -(d2*m1*pow(a1,2) - r[3]*sin(2*r[0] - 2*r[1])*pow(a2,2)*pow(l1,2)*pow(m2,2) + d2*cos(r[0] - r[1])*a2*l1*m2 + d2*pow(l1,2)*m2 + J1*d2)/(J1*J2 + pow(a2,2)*pow(l1,2)*pow(m2,2) + J2*pow(a1,2)*m1 + J1*pow(a2,2)*m2 + J2*pow(l1,2)*m2 - pow(a2,2)*pow(l1,2)*pow(m2,2)*pow(cos(r[0] - r[1]),2) + pow(a1,2)*pow(a2,2)*m1*m2);

        drdt[0] = L[0][0]*r[0] + L[0][1]*r[1] + L[0][2]*r[2] + L[0][3]*r[3];
        drdt[1] = L[1][0]*r[0] + L[1][1]*r[1] + L[1][2]*r[2] + L[1][3]*r[3];
        drdt[2] = L[2][0]*r[0] + L[2][1]*r[1] + L[2][2]*r[2] + L[2][3]*r[3];
        drdt[3] = L[3][0]*r[0] + L[3][1]*r[1] + L[3][2]*r[2] + L[3][3]*r[3];
    };
    solver(radODE, r, u);
};

auto pendulumAddO = [](SymbolicSet* O) -> void {

    ;
};

auto pendulumAddG = [](SymbolicSet* G) -> void {
    double V[4*4] = {0.247, 0.153, -0.023, -0.026,
                     0.153, 0.106, 0.026, -0.023,
                     -0.023, 0.026, 8.24, 3.893,
                     -0.026, -0.023, 3.893, 1.922};
    double c[4] = {M_PI, 0, 0, 0};
    G->addEllipsoid(V, c, INNER);
};

auto pendulumAddI = [](SymbolicSet* I) -> void {
    double q[4] = {M_PI, M_PI, 0, 0};
    I->addPoint(q);
};

int main() {

    double lbX[dimX]={M_PI / 2, 0, -5.7, -5.7};
    double ubX[dimX]={M_PI + 0.1, 2*M_PI, 5.7, 5.7};
    double lbU[dimU]={-3.5 * 9.81};
    double ubU[dimU]= {3.5 * 9.81};
    double etaU[dimU]= {17.1675};
    double etaRatio[dimX] = {3, 3, 3, 3};
    double tauRatio = 1;
    int nint = 5;

    double etaX[dimX]= {(M_PI/2+0.1)/118, 2*M_PI/118, 11.4/118, 11.4/118};
    double tau = 0.01;
    int numAbs = 1;

    int readXX = 0; // if X or U has changed, needs to be 0
    int readAbs = 0; // if X or U or O or dynamics has changed, needs to be 0

    Adaptive<X_type, U_type> abs(dimX, lbX, ubX, etaX, tau,
                                 dimU, lbU, ubU, etaU,
                                 etaRatio, tauRatio, nint,
                                 numAbs, readXX, readAbs);

    abs.initializeReach(pendulumAddG, pendulumAddI, pendulumAddO);
    abs.computeAbstractions(sysNext, radNext);

    int startAbs = 0;
    int minToGoCoarser = 3;
    int minToBeValid = 3;
//    abs.reach(startAbs, minToGoCoarser, minToBeValid, 1);
    abs.reachSCOTS();
}
