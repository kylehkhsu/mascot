#ifndef ADAPTIVE_HH_
#define ADAPTIVE_HH_

#include <string>
#include <vector>

#include "SymbolicModelGrowthBound.hh"
#include "FixedPoint.hh"
#include "RungeKutta4.hh"
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>


using std::cout;
using std::string;
using std::vector;
namespace scots {

template<class X_type, class U_type>
class Adaptive {
public:
    Cudd ddmgr_;

    // X and U
    int dimX_;
    double* lbX_;
    double* ubX_;
    double* etaX_;
    double tau_;
    int dimU_;
    double* lbU_;
    double* ubU_;
    double* etaU_;

    double etaRatio_ = 2;
    double tauRatio_ = 1.5;

    /* constructor */
    Adaptive(int dimX, double* lbX, double* ubX, double* etaX, double tau,
             int dimU, double* lbU, double* ubU, double* etaU) {
        dimX_ = dimX;
        lbX_ = lbX;
        ubX_ = ubX;
        etaX_ = etaX;
        tau_ = tau;
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
    void reach(sys_type sysNext, rad_type radNext, obst_type addO,
               goal_type addG, init_type addI, SymbolicSet* Fp = nullptr, int iter = 1) {

        SymbolicSet X(ddmgr_, dimX_, lbX_, ubX_, etaX_, tau_);
        SymbolicSet O(X);
        SymbolicSet G(X);
        SymbolicSet I(X);
        SymbolicSet U(ddmgr_, dimU_, lbU_, ubU_, etaU_, 0);
        SymbolicSet C(X,U);
        SymbolicSet X2(X, 1);

        X.addGridPoints();
        addO(&O);
        addG(&G);
        X.symbolicSet_ &= !O.symbolicSet_;

        if (Fp) {
            SymbolicSet Fn(X);
            innerApproximation(Fp, &Fn);
            G.symbolicSet_ |= Fn.symbolicSet_;
        }

        G.symbolicSet_ &= X.symbolicSet_;

        addI(&I);

        U.addGridPoints();

        cout << "\n\niteration: " << iter << '\n';
        printEta();
        cout << "tau: " << tau_ << '\n';

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

        if ((F.symbolicSet_ & I.symbolicSet_) != (ddmgr_.bddZero())) {
            X.writeToFile("X.bdd");
            O.writeToFile("O.bdd");
            I.writeToFile("I.bdd");
            return;
        }
        else {
            iter++;
            for (int i = 0; i < dimX_; i++) {
                etaX_[i] /= 2;
            }
            tau_ /= 1.5;
//            for (size_t i = 0; i < dimU_; i++) {
//                etaU_[i] /= 1.2;
//            }
            reach(sysNext, radNext, addO, addG, addI, &F, iter);
        }

    }

    /* function: test */
    template <class G_type>
    void test(G_type addG) {
        // X1 is coarser
        SymbolicSet X1(ddmgr_, dimX_, lbX_, ubX_, etaX_, tau_);
        for (int i = 0; i < dimX_; i++) {
            etaX_[i] /= etaRatio_;
        }
        tau_ /= tauRatio_;

        // X2 is finer
        SymbolicSet X2(ddmgr_, dimX_, lbX_, ubX_, etaX_, tau_);

        SymbolicSet X12(X1, X2);
        X1.addGridPoints();
        X2.addGridPoints();
        mapAbstractions(&X2, &X12);

        X1.printInfo(1);
        X2.printInfo(1);
        X12.printInfo(1);

        SymbolicSet F2(X2);
        addG(&F2);
        F2.printInfo(1);

        SymbolicSet F1(X1);
        innerCoarser(&F1, &F2, &X12);
        F1.printInfo(1);

        checkMakeDir("T");
        F1.writeToFile("T/F1.bdd");
        F2.writeToFile("T/F2.bdd");

    }

    /* function: innerCoarser */
    void innerCoarser(SymbolicSet* Fc, SymbolicSet* Ff, SymbolicSet* X12) {
        SymbolicSet Q12(*X12);
        Q12.symbolicSet_ = X12->symbolicSet_ & Ff->symbolicSet_;
        SymbolicSet R1(*Fc);
        R1.symbolicSet_ = Q12.symbolicSet_.ExistAbstract(Ff->getCube()); // & S1
        int* R1Minterm;
        SymbolicSet C12(*X12);
        BDD result = ddmgr_.bddZero();
        for (R1.begin(); !R1.done(); R1.next()) {
            R1Minterm = (int*)R1.currentMinterm();
            BDD R1BDD = R1.mintermToBDD(R1Minterm);
            C12.symbolicSet_ = Q12.symbolicSet_ & R1BDD;
            if (C12.symbolicSet_.CountMinterm(C12.nvars_) == dimX_ * dimX_) {
                result |= R1BDD;
            }
        }
        Fc->symbolicSet_ = result;
    }

    /* function: mapAbstractions */
    void mapAbstractions(SymbolicSet* Xf, SymbolicSet* X12) {

        int* XfMinterm;
        double Xpoint[dimX_] = {0};
        double XXpoint[dimX_ + dimX_] = {0};
        for (Xf->begin(); !Xf->done(); Xf->next()) {
            XfMinterm = (int*)Xf->currentMinterm();
            Xf->mintermToElement(XfMinterm, Xpoint);
            for (int i = 0; i < dimX_+dimX_; i++) {
                XXpoint[i] = Xpoint[i % dimX_];
            }
            X12->addPoint(XXpoint);
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

//    /* function: test */
//    void test() {
//        cout << lbX_[1] << '\n';
//    }

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

    /* function: checkMakeDir */
    template <class dir_type>
    void checkMakeDir(dir_type dirName) {
        DIR* dir = opendir(dirName);
        if (dir) {
            closedir(dir);
        }
        else if (ENOENT == errno) {
            int result = mkdir(dirName, 0777);
            (void) result;
        }
    }

};


}

#endif /* ADAPTIVE_HH_ */

