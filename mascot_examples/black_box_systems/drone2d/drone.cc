/*
 * drone.cc
 *
 *  created on: 07.11.2019
 *      author: kaushik
 */

/*
 * Quadcopter dynamics taken from the book https://www.kth.se/polopoly_fs/1.588039.1441112632!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf
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


#define dimX 2 /* X,Y coordinates */
#define dimX_hidden 10
#define dimU 4

/* data types for the ode solver */
typedef std::array<double, dimX> X_type;
typedef std::array<double, dimX+dimX_hidden> X_type_full;
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
auto  drone_post = overload
    (
        [](X_type &y, U_type &u, OdeSolver solver) -> void {
            /* The system is simulated only using the observed states.
             * The hidden ones are initialized randomly within some bound */
            X_type_full x;
            for (int i=0; i<dimX; i++)
                x[i]=y[i];
            /* bounds on hidden state variables (for initialization) */
            double lbX[dimX_hidden] = {0,0,0,0,-0.5,0,0,0,0,};
            double ubX[dimX_hidden] = {5.0,2*M_PI,2*M_PI,2*M_PI,1.0,1.0,1.0,0.1,0.1,0.1};
            /* randomly initialize the hidden state variables.
            * the hidden dimension variables are 3 to 11. */
            for (int i=0; i<dimX_hidden; i++) {
                x[dimX+i]=(rand() % 100*ceil(ubX[i]-lbX[i]))*0.01 + lbX[i];
//                x[dimX+i]=0;
            }
//            X_type x_curr = x;
            auto drone_ode = [](X_type_full &dxdt, X_type_full &x, const U_type &u) -> void {
                /* parameters */
                double m=1; /* in kg */
                double g=9.81; /* m/s */
                double b=54.2 * 1e-6; /* N/s/s */
                double l=0.24; /* m */
                double d=1.1*1e-6; /* Nm/s/s */
                double I_x=8.1*1e-3; /* Nms^2 */
                double I_y=8.1*1e-3; /* Nms^2 */
                double I_z=14.2*1e-3; /* Nms^2 */
                /* compute input forces and torques as function of the rotor speed (given by u) */
                double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
                double tau_x=b*l*(u[2]*u[2]-u[0]*u[0]);
                double tau_y=b*l*(u[3]*u[3]-u[1]*u[1]);
                double tau_z=d*(u[3]*u[3]+u[1]*u[1]-u[2]*u[2]-u[0]*u[0]);
                /* the angles are always between 0 and 2*M_PI */
                if (x[3]>2*M_PI) {
                  double rem = fmod(x[3],2*M_PI);
                  x[3] = rem;
                }
                if (x[3]<0) {
                  double rem = fmod(-x[3],2*M_PI);
                  x[3] = 2*M_PI - rem;
                }
                if (x[4]>2*M_PI) {
                  double rem = fmod(x[4],2*M_PI);
                  x[4] = rem;
                }
                if (x[4]<0) {
                  double rem = fmod(-x[4],2*M_PI);
                  x[4] = 2*M_PI - rem;
                }
                if (x[5]>2*M_PI) {
                  double rem = fmod(x[5],2*M_PI);
                  x[5] = rem;
                }
                if (x[5]<0) {
                  double rem = fmod(-x[5],2*M_PI);
                  x[5] = 2*M_PI - rem;
                }
                /* compute the derivative (eq. 2.25 in the book https://www.kth.se/polopoly_fs/1.588039.1441112632!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf) */
                dxdt[0] = x[6];
                dxdt[1] = x[7];
                dxdt[2] = x[8];
                dxdt[3] = x[10]*std::sin(x[5])/std::cos(x[4]) + x[11]*std::cos(x[5])/std::cos(x[4]);
                dxdt[4] = x[10]*std::cos(x[5]) - x[11]*std::sin(x[5]);
                dxdt[5] = x[9] + x[10]*std::sin(x[5])*std::tan(x[4]) + x[11]*std::cos(x[5])*std::tan(x[4]);
                dxdt[6] =-1/m *(std::sin(x[5])*std::sin(x[3])+std::cos(x[5])*std::cos(x[3])*std::sin(x[4])) * f_t;
                dxdt[7] =-1/m *(std::cos(x[3])*std::sin(x[5])-std::cos(x[5])*std::sin(x[3])*std::sin(x[4])) *f_t;
                dxdt[8] =-1/m *std::cos(x[5])*std::cos(x[4])*f_t+g;
                dxdt[9] =((I_y-I_z)/I_x)*x[10]*x[11] + 1/I_x*tau_x;
                dxdt[10]=((I_z-I_x)/I_y)*x[9]*x[11] + 1/I_y*tau_y;
                dxdt[11]=((I_x-I_y)/I_z)*x[9]*x[10] +1/I_z *tau_z;
            };
            solver(drone_ode, x, u);
            /* Send back the measured variables.
             * Measured variable dimensions are 0 to 2. */
            for (int i=0; i<dimX; i++) {
                y[i]=x[i];
            }
//            /* this is a hack to ignore the out-of-bound behaviors of all the variables other than the first 3 state variables (X, Y, Z coordinates) */
//        for (int i=3; i<dimX; i++) {
//            x[i] = x_curr[i];
//        }
            
        },
         [](X_type &y, U_type &u, OdeSolver solver, const std::vector<std::vector<double>> obs, std::vector<double>& unsafeAt) -> void {
        /* In this case the system is simulated using the full state space variables.
         * Initially, the hidden states are chosen randomly, and the next state is stored in a .txt file.
         * At each time step, the states are read from a .txt file, the system is simulated for one time step, and the
         * resulting states are appended to the same .txt file. */
        X_type_full x;
        std::ifstream file_read;
        file_read.open("system_trajectory.txt");
        /* flag set to true when the current state was read from file */
        bool state_read_from_file = false;
        if (file_read.is_open()) {
            std::string line;
            /* find the last position where x occurs in the file */
            int last_x_pos=-1;
            while(std::getline(file_read,line)) {
                if(line.find("x")!=std::string::npos) {
                    last_x_pos=file_read.tellg();
                }
            }
            /* the current state corresponds to the last occurrence of x */
            if (last_x_pos!=-1) {
                file_read.seekg(last_x_pos, ios::beg);
                /* read current state value in the next line and store in x */
                if(std::getline(file_read,line)) {
                    std::stringstream stream(line);
                    for (int j=0; j<dimX+dimX_hidden; j++)
                        stream >> x[j];
                    state_read_from_file=true;
                }
            }
        }
        file_read.close();
        /* in case the state couldn't be read from file, initialize the state and add to the file */
        if (!state_read_from_file) {
            /* no match found */
            /* initialize the measured dimesnion as per y */
            for (int i=0; i<dimX; i++)
                x[i]=y[i];
            /* bounds on hidden state variables (for initialization) */
            double lbX[dimX_hidden] = {0,0,0,0,-0.5,0,0,0,0,};
            double ubX[dimX_hidden] = {5.0,2*M_PI,2*M_PI,2*M_PI,1.0,1.0,1.0,0.1,0.1,0.1};
            /* randomly initialize the hidden state variables.
            * the hidden dimension variables are 3 to 11. */
            for (int i=dimX; i<dimX+dimX_hidden; i++) {
                x[i]=(rand() % 100*ceil(ubX[i]-lbX[i]))*0.01 + lbX[i];
            }
            /* add the initial state to the file */
            std::ofstream file_write;
            file_write.open("system_trajectory.txt",std::fstream::out);
            file_write << "x\n";
            for (int i=0; i<dimX+dimX_hidden; i++)
                file_write << x[i];
            file_write.close();
        }
         auto drone_ode = [](X_type_full &dxdt, X_type_full &x, const U_type &u) -> void {
             /* parameters */
             double m=1; /* in kg */
             double g=9.81; /* m/s */
             double b=54.2 * 1e-6; /* N/s/s */
             double l=0.24; /* m */
             double d=1.1*1e-6; /* Nm/s/s */
             double I_x=8.1*1e-3; /* Nms^2 */
             double I_y=8.1*1e-3; /* Nms^2 */
             double I_z=14.2*1e-3; /* Nms^2 */
             /* compute input forces and torques as function of the rotor speed (given by u) */
             double f_t=b*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
             double tau_x=b*l*(u[2]*u[2]-u[0]*u[0]);
             double tau_y=b*l*(u[3]*u[3]-u[1]*u[1]);
             double tau_z=d*(u[3]*u[3]+u[1]*u[1]-u[2]*u[2]-u[0]*u[0]);
             /* the angles are always between 0 and 2*M_PI */
             if (x[3]>2*M_PI) {
               double rem = fmod(x[3],2*M_PI);
               x[3] = rem;
             }
             if (x[3]<0) {
               double rem = fmod(-x[3],2*M_PI);
               x[3] = 2*M_PI - rem;
             }
             if (x[4]>2*M_PI) {
               double rem = fmod(x[4],2*M_PI);
               x[4] = rem;
             }
             if (x[4]<0) {
               double rem = fmod(-x[4],2*M_PI);
               x[4] = 2*M_PI - rem;
             }
             if (x[5]>2*M_PI) {
               double rem = fmod(x[5],2*M_PI);
               x[5] = rem;
             }
             if (x[5]<0) {
               double rem = fmod(-x[5],2*M_PI);
               x[5] = 2*M_PI - rem;
             }
             /* compute the derivative (eq. 2.25 in the book https://www.kth.se/polopoly_fs/1.588039.1441112632!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf) */
             dxdt[0] = x[6];
             dxdt[1] = x[7];
             dxdt[2] = x[8];
             dxdt[3] = x[10]*std::sin(x[5])/std::cos(x[4]) + x[11]*std::cos(x[5])/std::cos(x[4]);
             dxdt[4] = x[10]*std::cos(x[5]) - x[11]*std::sin(x[5]);
             dxdt[5] = x[9] + x[10]*std::sin(x[5])*std::tan(x[4]) + x[11]*std::cos(x[5])*std::tan(x[4]);
             dxdt[6] =-1/m *(std::sin(x[5])*std::sin(x[3])+std::cos(x[5])*std::cos(x[3])*std::sin(x[4])) * f_t;
             dxdt[7] =-1/m *(std::cos(x[3])*std::sin(x[5])-std::cos(x[5])*std::sin(x[3])*std::sin(x[4])) *f_t;
             dxdt[8] =-1/m *std::cos(x[5])*std::cos(x[4])*f_t+g;
             dxdt[9] =((I_y-I_z)/I_x)*x[10]*x[11] + 1/I_x*tau_x;
             dxdt[10]=((I_z-I_x)/I_y)*x[9]*x[11] + 1/I_y*tau_y;
             dxdt[11]=((I_x-I_y)/I_z)*x[9]*x[10] +1/I_z *tau_z;
         };
         solver(drone_ode, x, u, obs, unsafeAt);
        /* return the measured variables */
        for (int i=0; i<dimX; i++)
            y[i]=x[i];
        /* write the full state information in the file */
        std::ofstream file_write;
        file_write.open("system_trajectory.txt",std::fstream::app);
        file_write << "\n x\n";
        for (int i=0; i<dimX+dimX_hidden; i++)
            file_write << x[i];
        file_write.close();
         }
     );

/* because the system is treated as a black box, the growth bound is 0 in all dimensions */
auto radius_post = [](X_type &r, U_type &u, OdeSolver solver) -> void {
    r[0] = 0;
    r[1] = 0;
};


template<std::size_t SIZE_o>
auto spawnO = [](std::vector<std::array<double,SIZE_o*dimX>>& HO, std::vector<std::array<double,SIZE_o>>& ho, int verbose=0) -> void {
    /* pick a random set of obstacles in the X-Y axis region [2.4,2.5] X [0,5], modeled as number of thin boxes */
    int p = 20;
    for (int i=0; i<10; i++) {
        int toss = rand() % 100; /* generate a random number between 0 to 99 */
        if (toss<p) { /* with probability 0.01*p, each obstacle is actually present */
            /* intialize the H arrays */
            std::array<double,4*dimX> H={-1, 0,
                                            1, 0,
                                            0,-1,
                                            0, 1};
            std::array<double,4> box1 = {-2.5, 3.0, -0.5*i, 0.5*(i+1)};
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
    double side = 0.7; /* length of the side of the goal */
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

/* generate a single random initial state */
template<std::size_t SIZE_i>
auto generateInitial = [](std::vector<double>& init, const std::vector<std::array<double,SIZE_i>> hi) -> void {
    double toss1, toss2;
    toss1 = 0.01*(rand() % (int)(100*(hi[0][1]+hi[0][0])))-hi[0][0];
    toss2 = 0.01*(rand() % (int)(100*(hi[0][3]+hi[0][2])))-hi[0][2];
    init.push_back(toss1);
    init.push_back(toss2);
    for (int i=2; i< dimX; i++) {
        double x = eps;
        init.push_back(x);
    }
};


/****************************************************************************/
/* main computation */
/****************************************************************************/
int main() {

    double lbX[dimX] = {0, 0};
    double ubX[dimX] = {5.0, 5.0};
    
    double lbU[dimU] = {0,0,0,0};
    double ubU[dimU] = {300,300,300,300};
    double etaU[dimU] = {50,50,50,50};
    
    int nSubInt = 5;
    int systemNSubInt = 4;
    
    double etaX[dimX] = {0.1, 0.1};
    double tau = 0.1; /* must be an integer multiple of systemTau */
    double systemTau = 0.0005; /* time step for simulating the system trajectory */
    
    double etaRatio[dimX] = {2, 2};
    double tauRatio = 2; /* must be an integer */
    int p = 2;
    int verbose = 1;
    bool readTsFromFile = true;
    bool onlyRunTest = false;
    
    int numAbs = 3;
    
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
    srand(1573540980);
//    srand(1572389026);
    /* problematic seeds */
    // 1567743385, 1567744613, 1567750636(distance=-1 bug)
    
    
    /* number of samples used */
    int NN = 10;
    /* number of tests */
    int num_tests = 100;
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
    
    char logfile [] = "drone.log";
    double spec_final;
    TicToc timer;
    
    if (!onlyRunTest) {
        timer.tic();
        while (numAbs<=10) {
            spec_final = find_abst(x, u,
                                          drone_post, radius_post,
                                          lbX, ubX, lbU, ubU,
                                          etaX, etaU, tau, systemTau,
                                          numAbs, etaRatio, tauRatio,
                                          spawnO<4>, spawnG<4>, spawnI<4>,
                                          HO, HG, HI,
                                          ho, hg, hi, generateInitial<4>,
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
                     drone_post, radius_post,
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
