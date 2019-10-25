/*
 * vehicle.cc
 *
 *  created on: 01.09.2019
 *      author: kaushik
 */

/*
 * simple unperturbed vehicle
 *
 */

#include <array>
#include <iostream>
#include <cmath>
#include <time.h> /* time used to seed the random number */
#include <float.h> /* for the smallest positive number */
#define _USE_MATH_DEFINES

#include "BlackBoxReachWrapper.hh"

using namespace std;
using namespace scots;
using namespace helper;


#define dimX 2
#define dimU 2

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
            auto vehicle_ode = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
                dxdt[0] = u[0];
                dxdt[1] = u[1];
                //debug
        //        int r = rand()%10;
        //        dxdt[0] = r+u[0];
        //        dxdt[1] = r+u[1];
                // debug end
            };
            solver(vehicle_ode, x, u);
            // simple sampled time test trajectory
        //    x[0]+= u[0];
        //    x[1]+=u[1];
            
            /* saturate around the bounding box [0,5] X [0,5] */
//            double eps = 1e-13; /* very small number added to make sure that boundary cases are pessimistically resolved */
            if (x[0]<0)
                x[0]=0+eps;
            if (x[0]>5)
                x[0]=5-eps;
            
            if (x[1]<0)
                x[1]=0+eps;
            if (x[1]>5)
                x[1]=5-eps;
        },
         [](X_type &x, U_type &u, OdeSolver solver, const std::vector<std::vector<double>> obs, std::vector<double>& unsafeAt) -> void {
             auto vehicle_ode = [](X_type &dxdt, const X_type &x, const U_type &u) -> void {
                 dxdt[0] = u[0];
                 dxdt[1] = u[1];
                 //debug
                 //        int r = rand()%10;
                 //        dxdt[0] = r+u[0];
                 //        dxdt[1] = r+u[1];
                 // debug end
             };
             solver(vehicle_ode, x, u, obs, unsafeAt);
             // simple sampled time test trajectory
             //    x[0]+= u[0];
             //    x[1]+=u[1];
             
             /* saturate around the bounding box [0,5] X [0,5] */
//             double eps = 1e-13; /* very small number added to make sure that boundary cases are pessimistically resolved */
             if (x[0]<0)
                 x[0]=0+eps;
             if (x[0]>5)
                 x[0]=5-eps;
             
             if (x[1]<0)
                 x[1]=0+eps;
             if (x[1]>5)
                 x[1]=5-eps;
             
             if (unsafeAt.size()!=0) {
                 if (unsafeAt[0]<0)
                     unsafeAt[0]=0+eps;
                 if (unsafeAt[0]>5)
                     unsafeAt[0]=5-eps;
                 
                 if (unsafeAt[1]<0)
                     unsafeAt[1]=0+eps;
                 if (unsafeAt[1]>5)
                     unsafeAt[1]=5-eps;
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
    /* pick a random set of obstacles in the region [2.4,2.5] X [0,5], modeled as number of thin boxes */
    int p = 30;
    for (int i=0; i<10; i++) {
        int toss = rand() % 100; /* generate a random number between 0 to 99 */
        if (toss<p) { /* with probability 0.01*p, each obstacle is actually present */
            /* intialize the H arrays */
            std::array<double,4*dimX> H={-1, 0,
                1, 0,
                0,-1,
                0, 1};
            std::array<double,4> box1 = {-2.4, 2.5, -0.5*i, 0.5*(i+1)};
            ho.push_back(box1);
            HO.push_back(H);
            if (verbose>0)
                cout << "Obstacle set: [" << -box1[0] << "," << box1[1] << "] X [" << -box1[2] << "," << box1[3] << "].\n";
        }
    }
};

template<std::size_t SIZE_g>
auto spawnG = [](std::vector<std::array<double,SIZE_g*dimX>>& HG, std::vector<std::array<double,SIZE_g>>& hg, int verbose=0) -> void {
    double toss1, toss2, toss3;
    /* pick a random goal set */
    /* the single goal set is a square of fixed size placed anywhere in the band [0,2] X [0,5] U [3,5] X [0,5] with some probability distribution */
    /* intialize the H arrays */
    std::array<double,4*dimX> H={-1, 0,
        1, 0,
        0,-1,
        0, 1};
    double side = 1.95; /* length of the side of the goal */
    toss1 = rand() % 100; /* generate a random number between 0 to 99 */
    if (toss1<98) { /* with 98% probability, it is in [3,5] X [0,5] */
        /* randomly generate the lower left corner of the square goal */
        toss2 = (rand() % (20-(int)side*10))*0.1 + 3; /* generate a random fraction between 3 to (5-side) */
    } else {
        /* randomly generate the lower left corner of the square goal */
        toss2 = (rand() % (20-(int)side*10))*0.1; /* generate a random fraction between 0 to (2-side) */
    }
    toss3 = (rand() % (50-(int)side*10))*0.1; /* generate a random fraction between 0 to (5-side) */
    std::array<double,4> box2 = {-toss2, toss2+side, -toss3, toss3+side};
    hg.push_back(box2);
    HG.push_back(H);
    if (verbose>0)
        cout << "Goal set: [" << -box2[0] << "," << box2[1] << "] X [" << -box2[2] << "," << box2[3] << "].\n";
};

template<std::size_t SIZE_i>
auto spawnI = [](std::vector<std::array<double,4*dimX>>& HI, std::vector<std::array<double,SIZE_i>>& hi,int verbose=0) -> void {
    double toss1, toss2, toss3;
    /* intialize the H arrays */
    std::array<double,4*dimX> H={-1, 0,
        1, 0,
        0,-1,
        0, 1};
    /* pick a random initial set */
    /* the initial state set is a small square of length side whose lower left corner is a random point in the region [0,2-side) X [0,5-side] U (3,5-side] X [0,5-side] */
    double side = 0.5;
    toss1 = rand() % 100; /* generate a random number between 0 to 99 */
    if (toss1<98) { /* with 98% probability, it is in [0,2-side] X [0,5-side] */
        /* randomly generate the x0 coordinate of the initial state */
        toss2 = (rand() % (int)(20-ceil(side*10)))*0.1;
    } else {
        /* randomly generate the x0 coordinate of the initial state */
        toss2 = (rand() % (int)(20-ceil(side*10)) + 1)*0.1 + 3;
    }
    toss3 = (rand() % (int)(50-ceil(side*10)))*0.1; /* generate a random fraction between 0 to 5-side */
    std::array<double,4> box3 = {-toss2, toss2+side, -toss3, toss3+side};
    hi.push_back(box3);
    HI.push_back(H);
    if (verbose>0)
        cout << "Initial set: [" << -box3[0] << "," << box3[1] << "] X [" << -box3[2] << "," << box3[3] << "].\n";
};


/****************************************************************************/
/* main computation */
/****************************************************************************/
int main() {

    double lbX[dimX] = {0, 0};
    double ubX[dimX] = {5.0, 5.0};
    
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
    int verbose = 1;
    bool readTsFromFile = false;
    
    int numAbs = 4;
    
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
    srand(1568788662);
//    srand(1568250847);
    /* problematic seeds */
    // 1567743385, 1567744613, 1567750636(distance=-1 bug)
    
    
    /* number of samples used */
    int NN = 100;
    /* number of tests */
    int num_tests = 50;
    /* allowed distance between an abstract trajectory and the unsafe states for the controllers good for the abstraction but bad for the system */
//    double allowedDistance = 0.15;
    X_type explRadius = {0.25,0.25};
    /* exploration horizon in time units */
    double explHorizon = 1;
    double reqd_success_rate[2] = {0.02, 0.9}; /* when to stop: when the fraction of refinement to the total iteration in loop 1 is below reqd_success_rate[0], and the success rate in loop 2 is above reqd_success_rate[1]. */
//    double reqd_success_rate = 0.9;
    
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
