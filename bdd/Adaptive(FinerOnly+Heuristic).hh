#ifndef ADAPTIVE_HH_
#define ADAPTIVE_HH_

#include <string>
#include <vector>

#include "SymbolicModelGrowthBound.hh"
#include "FixedPoint.hh"
#include "RungeKutta4.hh"


using std::cout;
using std::string;
using std::vector;
namespace scots {

template<class X_type, class U_type>
class Adaptive {
public:

    // X and U
    int dimX_;
    double* lbX_;
    double* ubX_;
    double* etaX_;
    int dimU_;
    double* lbU_;
    double* ubU_;
    double* etaU_;

    /* constructor */
    Adaptive(int dimX, double* lbX, double* ubX, double* etaX,
             int dimU, double* lbU, double* ubU, double* etaU) {
        dimX_ = dimX;
        lbX_ = lbX;
        ubX_ = ubX;
        etaX_ = etaX;
        dimU_ = dimU;
        lbU_ = lbU;
        ubU_ = ubU;
        etaU_ = etaU;
    }

    ~Adaptive() {
        ;
    }

    /* function: reach */
    template <class sys_type, class rad_type, class obst_type, class goal_type, class init_type>
    void reach(sys_type sysNext, rad_type radNext, obst_type addObstacles,
               goal_type constructGoal, init_type constructInit, SymbolicSet* Fp = nullptr, int iter = 1) {
        Cudd ddmgr;
        SymbolicSet X(ddmgr, dimX_, lbX_, ubX_, etaX_);
        X.addGridPoints();
        addObstacles(&X);

//        cout << "X:\n";
//        X.printInfo(1);

        SymbolicSet G(X);
        constructGoal(&G);
        // G.symbolicSet_ &= Fp
        if (Fp) {
            SymbolicSet Fn(X);
            innerApproximation(Fp, &Fn);
            G.symbolicSet_ |= Fn.symbolicSet_;
        }

        G.symbolicSet_ &= X.symbolicSet_;
//        cout << "G:\n";
//        G.printInfo(1);

        SymbolicSet I(X);
        constructInit(&I);

        SymbolicSet U(ddmgr, dimU_, lbU_, ubU_, etaU_);
        U.addGridPoints();

        SymbolicSet C(X, U);
        SymbolicSet X2(X, 1);

        cout << "\n\niteration: " << iter << '\n';
        printEta();

        SymbolicSet Xbar(X);
        Xbar.symbolicSet_ = ((!G.symbolicSet_) & X.symbolicSet_);

        SymbolicModelGrowthBound<X_type, U_type> abstraction(&Xbar, &U, &X2);
        abstraction.computeTransitionRelation(sysNext, radNext, X.symbolicSet_);
        std::cout << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;
        FixedPoint fp(&abstraction);
        C.symbolicSet_ = fp.reach(G.symbolicSet_, I.symbolicSet_, 1);

        string cStr = "C";
        cStr += std::to_string(iter);
        cStr += ".bdd";
        char cChar[20];
        size_t cLength = cStr.copy(cChar, cStr.length() + 1);
        cChar[cLength] = '\0';
        C.writeToFile(cChar);
//        cout << "C:\n";
//        C.printInfo(1);

        string gStr = "G";
        gStr += std::to_string(iter);
        gStr += ".bdd";
        char gChar[20];
        size_t gLength = gStr.copy(gChar, gStr.length() + 1);
        gChar[gLength] = '\0';
        G.writeToFile(gChar);

        SymbolicSet F(X);
        F.symbolicSet_ = C.symbolicSet_.ExistAbstract(U.getCube());

//        SymbolicSet T = abstraction.getTransitionRelation();
//        cout << "T:\n";
//        T.printInfo(1);

        if ((F.symbolicSet_ & I.symbolicSet_) != (ddmgr.bddZero())) {
            X.writeToFile("X.bdd");
            X.complement();
            X.writeToFile("X_obst.bdd");
            return;
        }
        else {
            iter++;
            for (int i = 0; i < dimX_; i++) {
                etaX_[i] /= 1.5;
            }
//            for (size_t i = 0; i < dimU_; i++) {
//                etaU_[i] /= 1.2;
//            }
            reach(sysNext, radNext, addObstacles, constructGoal, constructInit, &F, iter);
        }

    }

    /* function: preProcess */
    template <class sys_type, class rad_type, class no_type, class obst_type, class goal_type, class init_type>
    void preProcess(sys_type sysNext, rad_type radNext, no_type radNo, obst_type addObstacles, goal_type constructGoal, init_type constructInit) {
        Cudd ddmgr;
        SymbolicSet X(ddmgr, dimX_, lbX_, ubX_, etaX_);
        X.addGridPoints();
        addObstacles(&X);

        SymbolicSet G(X);
        constructGoal(&G);
        G.symbolicSet_ &= X.symbolicSet_;

        SymbolicSet I(X);
        constructInit(&I);

        SymbolicSet U(ddmgr, dimU_, lbU_, ubU_, etaU_);
        U.addGridPoints();

        SymbolicSet C(X, U);
        SymbolicSet X2(X, 1);

        SymbolicModelGrowthBound<X_type, U_type> abstractionPre(&X, &U, &X2);
        abstractionPre.computeTransitionRelation(sysNext, radNo, X.symbolicSet_);

        FixedPoint fpPre(&abstractionPre);
        C.symbolicSet_ = fpPre.reach(G.symbolicSet_, I.symbolicSet_, 1);

        SymbolicSet F(X);
        F.symbolicSet_ = C.symbolicSet_.ExistAbstract(U.getCube());
        if ((F.symbolicSet_ & I.symbolicSet_) != (ddmgr.bddZero())) {

            C.writeToFile("Cpre.bdd");


            SymbolicSet Xheur(X);
            Xheur.printInfo(1);

            X_type xPoint;
            U_type uPoint;

            I.begin();
            int* minterm = (int*) I.currentMinterm();
            I.mintermToElement(minterm, &xPoint[0]);

            vector<double> xVector;
            vector<size_t> ind;
            double H[dimX_ * dimX_] = {0.25, 0,
                                       0, 0.25};
            double c[dimX_] = {0};

            for (size_t i = 0; i < dimX_; i++) {
                xVector.push_back(xPoint[i]);
                ind.push_back(i);
                c[i] = xPoint[i];
            }
            Xheur.addEllipsoid(H, c, INNER);

            while (!(G.isElement(xVector))) {
                vector<vector<double>> inputs = C.setValuedMap(xVector,ind);
//                cout << inputs.size() << '\n';
                vector<double> uVector = inputs[0];
                for (size_t i = 0; i < dimU_; i++) {
                    uPoint[i] = uVector[i];
                }
                printx(xPoint);
                printu(uPoint);

                sysNext(xPoint, uPoint);

                for (size_t i = 0; i < dimX_; i++) {
                    xVector[i] = xPoint[i];
                    c[i] = xPoint[i];
                }
                Xheur.addEllipsoid(H, c, INNER);

            }
            Xheur.printInfo(1);
            addObstacles(&Xheur);
            Xheur.printInfo(1);
            Xheur.writeToFile("Xheur.bdd");

            SymbolicModelGrowthBound<X_type, U_type> abstraction(&Xheur, &U, &X2);
            abstraction.computeTransitionRelation(sysNext, radNext, X.symbolicSet_);
            FixedPoint fp(&abstraction);
            C.symbolicSet_ = fp.reach(G.symbolicSet_, I.symbolicSet_, 1);
            F.symbolicSet_ = C.symbolicSet_.ExistAbstract(U.getCube());
            if ((F.symbolicSet_ & I.symbolicSet_) != (ddmgr.bddZero())) {
                C.writeToFile("C1.bdd");
                X.writeToFile("X.bdd");
                X.complement();
                X.writeToFile("X_obst.bdd");
                G.writeToFile("G1.bdd");
            }
            else {
                cout << "Heuristic tube failed to yield controller.\n";
            }

        }
        else {
            std::ostringstream os;
            os << "preProcess failed to generate heuristic path.\n";
            throw std::runtime_error(os.str().c_str());
        }


    }

    /* function: innerApproximation */
    void innerApproximation(SymbolicSet* Fp, SymbolicSet* Fn) {
        int dimX = dimX_;
        double* etaX = etaX_;
        auto f = [Fp, dimX, etaX](double* x)->bool {
            std::vector<double> exPlus(x, x + dimX);
            std::vector<double> exMinus(x, x + dimX);
            for (int i = 0; i < dimX; i++) {
                exPlus[i] += etaX[i];
                exMinus[i] -= etaX[i];
            }
            return (Fp->isElement(exPlus) && Fp->isElement(exMinus));
        };
        Fn->addByFunction(f);
    }

    /* function: printEta */
    void printEta() {
        cout << "X eta: ";
        printArray(etaX_, dimX_);
        cout << "U eta: ";
        printArray(etaU_, dimU_);
    }


    /* function: saveSystem */
    void saveSystem(SymbolicSet* X, SymbolicSet* G) {
        X->writeToFile("X.bdd");
        G->writeToFile("G.bdd");
    }

    /* function: test */
    void test() {
        cout << lbX_[1] << '\n';
    }

    /* function: printx */
    void printx (X_type x) {
        cout << "x: ";
        for (size_t i = 0; i < x.size(); i++) {
            cout << x[i] << ' ';
        }
        cout << '\n';
    }

    /* function: printu */
    void printu (U_type u) {
        cout << "u: ";
        for (size_t i = 0; i < u.size(); i++) {
            cout << u[i] << ' ';
        }
        cout << '\n';
    }

    /* function: printArray
            * prints std::array or int array or double array */
    template <class T>
    void printArray (T array) {
        for (size_t i = 0; i < array.size(); i++) {
            cout << array[i] << ' ';
        }
        cout << '\n';
    }
    void printArray (int* array, size_t size) {
        for (size_t i = 0; i < size; i++) {
            cout << array[i] << ' ';
        }
        cout << '\n';
    }
    void printArray (double* array, size_t size) {
        for (size_t i = 0; i < size; i++) {
            cout << array[i] << ' ';
        }
        cout << '\n';
    }

};


}

#endif /* ADAPTIVE_HH_ */
