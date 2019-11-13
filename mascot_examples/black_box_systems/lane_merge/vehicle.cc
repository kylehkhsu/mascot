/*
 * vehicle.cc
 *
 *  created on: 15.09.2019
 *      author: kaushik
 */

/*
 * controller synthesis for lane merging
 *
 */

#include <array>
#include <iostream>
#include <cmath>
#include <time.h> /* time used to seed the random number */
#include <float.h> /* for the smallest positive number */

/* headers for chrono */
#include <ChDriver.h>
#include <ChVehicle.h>
#include <Sedan_Vehicle.h>
#include <ChCoordsys.h>
#include <ChVector.h>

#define _USE_MATH_DEFINES

#include "BlackBoxReachWrapper.hh"

//using namespace std;
using namespace scots;
using namespace helper;
using namespace chrono;


#define dimX 5 // for 1 other car in the first lane
#define dimU 2

#define lane_width 3.75 // lane width in metres
#define accl_lane_len 200 // length of acceleration lane in metres
#define min_lateral_clearance 2
#define min_longitudinal_clearance 40
#define min_prescribed_speed 80
#define max_speed 130

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimU> U_type;

/* overloading the system post function with 3 and 5 arguments. this trick is taken from here: http://martinecker.com/martincodes/lambda-expression-overloading/ */
template <class F1, class F2>
struct overload_set : F1, F2
{
    overload_set(F1 f1, F2 f2)
    : F1(f1), F2(f2)
    {}
    
    using F1::operator();
    using F2::operator();
};

template <class F1, class F2>
overload_set<F1, F2> overload(F1 f1, F2 f2)
{
    return overload_set<F1, F2>(f1, f2);
}

/* we consider a sampled time unperturbed system  */
auto  vehicle_post = overload
    (
        [](X_type &x, U_type &u, OdeSolver solver) -> void {
            /* normal vehicle post */
            ChSystem ch_sys();
            vehicle::sedan::Sedan_Vehicle sedan(&ch_sys);
            ChVector<Real> pos(x[0],x[1],0.0);
            sedan.Initialize(pos);
            sedan.Synchronize();
            sedan.setStepsize(solver.getStepSize());
            sedan.Advance(solver.tau_);
            ChVector<Real> pos2 = sedan.getVehiclePos();
            x[0] = pos2.x();
            x[1] = pos2.y();
            x[2] = sedan.getVehicleSpeed();
            
            /* other vehicle post (assuming constant velocity) */
            x[3] += x[4]*solver.tau_;
        },
         [](X_type &x, U_type &u, OdeSolver solver, const std::vector<std::vector<double>> obs, std::vector<double>& unsafeAt) -> void {
                //vehicle post with collision check
             }
         }
     );

/* because the system is treated as a black box, the growth bound is 0 in all dimensions */
auto radius_post = [](X_type &r, U_type &u, OdeSolver solver) -> void {
    r[0] = 0;
    r[1] = 0;
};


template<std::size_t SIZE_o>
auto spawnO = [](std::vector<std::array<double,SIZE_o*dimX>>& HO, std::vector<std::array<double,SIZE_o>>& ho, int verbose=0) -> void {
    /* fixed obstacle: end of lane */
    std::array<double,4*dimX> H1={-1, 0, 0, 0, 0,
                                   1, 0, 0, 0, 0,
                                   0, -1, 0, 0, 0,
                                   0, 1, 0, 0, 0};
    std::array<double,4> h1 = {-lane_width, 2*lane_width, 0, accl_lane_len};
    ho.push_back(h1);
    HO.push_back(H1);
    /* distance from the other car */
    /* using the polytope: (x1-y1) > 2 & (y1-x1) >2 & (x2-y2) > 40 & (y2-x2) > 40 */
    std::array<double,4*dimX> H2={-1, 0, 0, 1, 0,
        1, 0, 0, -1, 0,
        0, -1, 0, 0, 1,
        0, 1, 0, 0, -1};
    std::array<double,4> h2 = {-min_lateral_clearance, -min_lateral_clearance, -min_longitudinal_clearance, -min_longitudinal_clearance};
    ho.push_back(h2);
    HO.push_back(H2);
    /* minimum velocity in the first lane */
    std::array<double,4*dimX> H3={-1, 0, 0, 0, 0,
        1, 0, 0, 0, 0,
        0, 0, -1, 0, 0,
        0, 0, 1, 0, 0};
    std::array<double,4> h3 = {0, lane_width, 0, min_prescribed_speed};
    ho.push_back(h3);
    HO.push_back(H3);
};

template<std::size_t SIZE_g>
auto spawnG = [](std::vector<std::array<double,SIZE_g*dimX>>& HG, std::vector<std::array<double,SIZE_g>>& hg, int verbose=0) -> void {
    /* to reach the 1st lane before the end of acceleration lane with a certain minimum speed */
    std::array<double,5*dimX> H={-1, 0, 0, 0, 0,
        1, 0, 0, 0, 0,
        0, -1, 0, 0, 0,
        0, 0, -1, 0, 0,
        0, 0, 1, 0, 0};
    std::array<double,5> h = {0, lane_width, -accl_lane_len, -min_prescribed_speed, max_speed};
    ho.push_back(h);
    HO.push_back(H);
};

template<std::size_t SIZE_i>
auto spawnI = [](std::vector<std::array<double,4*dimX>>& HI, std::vector<std::array<double,SIZE_i>>& hi,int verbose=0) -> void {
    double self_speed = (rand() % 100)*0.1 + 50;
    double other_car_pos = (rand() % (std::ceil(accl_lane_len*10)))*0.1;
    double other_car_speed = (rand() % (std::ceil(max_speed-min_prescribed_speed)*10))*0.1 + min_prescribed_speed;
    std::array<double,10*dimX> H={-1, 0, 0, 0, 0,
        1, 0, 0, 0, 0,
        0, -1, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, -1, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, -1, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, -1,
        0, 0, 0, 0, 1};
    std::array<double,10> h = {-1.5*lane_width, 1.5*lane_width, 0, 0, -self_speed, self_speed, -other_car_pos, other_car_pos, -other_car_speed, other_car_speed};
    ho.push_back(h);
    HO.push_back(H);
};


/****************************************************************************/
/* main computation */
/****************************************************************************/
int main() {

    double lbX[dimX] = {0, 0, min_prescribed_speed-5, 0, min_prescribed_speed-5};
    double ubX[dimX] = {2*lane_width, accl_lane_len+10, max_speed, accl_lane_len+10, max_speed};
    
    double lbU[dimU] = {-2.0, -2.0};
    double ubU[dimU] = {2.0, 2.0};
    double etaU[dimU] = {0.2, 0.2};
    
    int nSubInt = 5;
    int systemNSubInt = 10;
    
    double etaX[dimX] = {0.3, 0.3};
    double tau = 0.25; /* must be an integer multiple of systemTau */
    double systemTau = 0.0005; /* time step for simulating the system trajectory */
    
    double etaRatio[dimX] = {2, 2};
    double tauRatio = 2; /* must be an integer */
    int p = 2;
    int verbose = 0;
    bool readTsFromFile = false;
    
    int numAbs = 6;
    
    /* setup colors for printing result */
    bool useColors = true; /* colored output not compatible with verbose > 0 so far */
//    const unsigned char colors[] =
//    {
//        0x02, /* foreground green */
//        0x04  /* foreground red */
//    };
//
//    HANDLE hstdout = GetStdHandle( STD_OUTPUT_HANDLE );
//    // Remember how things were when we started
//    CONSOLE_SCREEN_BUFFER_INFO csbi;
//    GetConsoleScreenBufferInfo( hstdout, &csbi );
    
    /* seed the random number generator with the current calendar value */
    int seed = time(NULL);
    cout << "\nSeed used for the random number generator : " << seed << "\n\n";
    srand(seed);
//    srand(1568531735);
//    srand(1568250847);
    /* problematic seeds */
    // 1567743385, 1567744613, 1567750636(distance=-1 bug)
    
    
    /* number of samples used */
    int NN = 1000;
    /* number of tests */
    int num_tests = 100;
    /* allowed distance between an abstract trajectory and the unsafe states for the controllers good for the abstraction but bad for the system */
//    double allowedDistance = 0.15;
    X_type explRadius = {0.25,0.25};
    /* exploration horizon in time units */
    double explHorizon = 1;
//    double reqd_success_rate[2] = {1, 0.9}; /* when to stop: when the fraction of unsound controller (existence but failure on system) in loop 1 is below reqd_success_rate[0], and the success rate in loop 2 is above reqd_success_rate[1]. */
    double reqd_success_rate = 0.9;
    
    /* maximum value of spec based on knowledge of the size of the goal
     * can be set to be the maximum side length of the state space if goal is not known precisely */
    double spec_max = 0.8;
    
    X_type x;
    U_type u;
    
    /* intialize the H arrays */
    std::vector<std::array<double,4*dimX>> HO;
    std::vector<std::array<double,4*dimX>> HG;
    std::vector<std::array<double,4*dimX>> HI;
    /* initialize the vectors containing the obstacles, goals, and the intial sets */
    std::vector<std::array<double,4>> ho;
    std::vector<std::array<double,4>> hg;
    std::vector<std::array<double,4>> hi;
    
    char logfile [] = "vehicle.log";
    double spec_final;
    TicToc timer;
    
    if (!readTsFromFile) {
        timer.tic();
        while (numAbs<=10) {
            spec_final = find_abst(x, u,
                                          vehicle_post, radius_post,
                                          lbX, ubX, lbU, ubU,
                                          etaX, etaU, tau, systemTau,
                                          numAbs, etaRatio, tauRatio,
                                          spawnO<4>, spawnG<4>, spawnI<4>,
                                          HO, HG, HI,
                                          ho, hg, hi,
                                          nSubInt, systemNSubInt, p,
                                          NN, explRadius, explHorizon, reqd_success_rate, spec_max,
                                          readTsFromFile, useColors, logfile, verbose);
            /* check if the number of abstractions need to be increased due to too high SPEC */
            if (spec_final==-1)
                numAbs++;
            else
                break;
        }

        cout << "\nFinal SPEC: " << spec_final << ".\n";


        cout << "\nTime taken by abstraction computation: " << timer.toc()/60 << " minutes\n";
    } else {
        cout << "\nEnter spec : ";
        cin >> spec_final;
    }
    
    /* use the computed abstraction for an actual controller synthesis */
    timer.tic();
    System sys(dimX, lbX, ubX, etaX, tau,
               dimU, lbU, ubU, etaU,
               etaRatio, tauRatio, nSubInt, numAbs);
    
    BlackBoxReach* abs= new BlackBoxReach(logfile,verbose);
    abs->initialize(&sys,systemTau,systemNSubInt);
    abs->loadTs();
    
    test_abstraction(abs, spec_final,
                     x, u,
                     vehicle_post, radius_post,
                     spawnO<4>, spawnG<4>, spawnI<4>,
                     HO, HG, HI,
                     ho, hg, hi,
                     p,
                     NN, num_tests,
                     useColors, logfile, verbose);
    
    cout << "\nTime taken by testing: " << timer.toc()/60 << " minutes.\n";
//    cout << "\nTime taken by the abstraction computation = " << comp_time << "s";
//    cout << "\nTime taken by the "
    
    
    return 1;
}
