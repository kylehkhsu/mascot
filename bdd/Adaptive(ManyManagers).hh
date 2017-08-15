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

    int dimX_;
    double* lbX_;
    double* ubX_;
    int dimU_;
    double* lbU_;
    double* ubU_;
    double* etaU_;
    double etaRatio_;
    double tauRatio_;
    int numAbs_;
    int startAbs_;

    vector<double*> etaX_;
    vector<double*> tau_;

    vector<Cudd*> ddmgrs_;
    vector<SymbolicSet*> Xs_;
    vector<SymbolicSet*> Us_;
    vector<SymbolicSet*> X2s_;

    vector<Cudd*> ddmgrsXX_;
    vector<SymbolicSet*> XXs_;




    /* constructor */
    template<class O_type, class sys_type, class rad_type>
    Adaptive(int dimX, double* lbX, double* ubX, double* etaX, double tau,
             int dimU, double* lbU, double* ubU, double* etaU,
             double etaRatio, double tauRatio,
             int numAbs, int startAbs,
             O_type addO, sys_type sysNext, rad_type radNext) {
        dimX_ = dimX;
        lbX_ = lbX;
        ubX_ = ubX;
        dimU_ = dimU;
        lbU_ = lbU;
        ubU_ = ubU;
        etaU_ = etaU;

        etaRatio_ = etaRatio;
        tauRatio_ = tauRatio;

        numAbs_ = numAbs;
        startAbs_ = startAbs;

        initializeEtaTau(etaX, tau);
        initializeDdmgrs();
        initializeXs();
        initializeUs();
        initializeX2s();

        initializeDdmgrsXX();
        initializeXXs();


        debug();
        test();

    }

    ~Adaptive() {
        deleteVecArray(etaX_);
        deleteVec(tau_);
        deleteVec(Xs_);
        deleteVec(Us_);
        deleteVec(X2s_);
        deleteVec(ddmgrs_); // has to be freed after SymbolicSets

        deleteVec(XXs_);
        deleteVec(ddmgrsXX_);

    }

    /* function: mapAbstractions */
    void mapAbstractions(SymbolicSet* Xf, SymbolicSet* Xcf) {
        int* XfMinterm;
        double Xpoint[dimX_] = {0};
        double XXpoint[dimX_ * 2] = {0};
        for (Xf->begin(); !Xf->done(); Xf->next()) {
            XfMinterm = (int*)Xf->currentMinterm();
            Xf->mintermToElement(XfMinterm, Xpoint);
            for (int i = 0; i < dimX_ * 2; i++) {
                XXpoint[i] = Xpoint[i % dimX_];
            }
            Xcf->addPoint(XXpoint);
        }
    }

    /* function: debug */
    void debug() {
        cout << "\n--------------------debug--------------------\n";
        printEtaX();
        printTau();
    }

    /* function: test */
    void test() {
        cout << "\n--------------------test--------------------\n";
        Cudd ddmgr1;
        Cudd ddmgr2;
        SymbolicSet X1(ddmgr1, dimX_, lbX_, ubX_, etaX_[0], 0);
        SymbolicSet X2(ddmgr2, dimX_, lbX_, ubX_, etaX_[0], 0);
        X1.addGridPoints();
        X1.printInfo(1);
        X2.symbolicSet_ = X1.symbolicSet_;
        X2.printInfo(1);

    }

    /* function: initializeEtaTau */
    void initializeEtaTau(double* etaX, double tau) {
        double etaCur[dimX_] = {0};
        for (int i = 0; i < dimX_; i++) {
            etaCur[i] = etaX[i];
        }
        double tauCur = tau;

        for (int i = 0; i < numAbs_; i++) {
            double* etai = new double[dimX_];
            double* taui = new double;
            for (int j = 0; j < dimX_; j++) {
                etai[j] = etaCur[j];
            }
            *taui = tauCur;
            etaX_.push_back(etai);
            tau_.push_back(taui);

            for (int j = 0; j < dimX_; j++) {
                etaCur[j] /= etaRatio_;
            }
            tauCur /= tauRatio_;
        }
    }


    /* function: initializeDdmgrs */
    void initializeDdmgrs() {
        for (int i = 0; i < numAbs_; i++) {
            Cudd* ddmgr = new Cudd;
            ddmgrs_.push_back(ddmgr);
        }
    }

    /* function: initializeXs */
    void initializeXs() {
        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* X = new SymbolicSet(*ddmgrs_[i], dimX_, lbX_, ubX_, etaX_[i], tau_[i][0]);
            X->addGridPoints();
            Xs_.push_back(X);
            X->printInfo(1);
        }
    }

    /* function: initializeUs */
    void initializeUs() {
        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* U = new SymbolicSet(*ddmgrs_[i], dimU_, lbU_, ubU_, etaU_, tau_[i][0]);
            U->addGridPoints();
            Us_.push_back(U);
            U->printInfo(1);
        }
    }

    /* function: initialize X2s */
    void initializeX2s() {
        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* X2 = new SymbolicSet(*Xs_[i], 1);
            X2->addGridPoints();
            X2s_.push_back(X2);
            X2->printInfo(1);
        }
    }

    /* function: initializeDdmgrsXX */
    void initializeDdmgrsXX() {
        for (int i = 0; i < numAbs_ - 1; i++) {
            Cudd* ddmgrXX = new Cudd;
            ddmgrsXX_.push_back(ddmgrXX);
        }
    }

    /* function: initializeXXs */
    void initializeXXs() {
//        for (int i = 0; i < numAbs_ - 1; i++) {
//            SymbolicSet Xc(*ddmgrsXX_[i], dimX_, lbX_, ubX_, etaX_[i], tau_[i][0]);
//            SymbolicSet Xf(*ddmgrsXX_[i], dimX_, )
//        }
        ;
    }


    /* function: deleteVec */
    template<class vec_type>
    void deleteVec(vec_type vec) {
        for (size_t i = 0; i < vec.size(); i++) {
            delete vec[i];
        }
    }

    /* function: deleteVecArray */
    template<class vec_type>
    void deleteVecArray(vec_type vec) {
        for (size_t i = 0; i < vec.size(); i++) {
            delete[] vec[i];
        }
    }

    /* function: printEtaX */
    void printEtaX() {
        cout << "etaX_:\n";
        for (size_t i = 0; i < etaX_.size(); i++) {
            cout << "abstraction " << i << ": ";
            for (int j = 0; j < dimX_; j++) {
                cout << etaX_[i][j] << " ";
            }
            cout << '\n';
        }
    }

    /* function: printTau */
    void printTau() {
        cout << "tau_:\n";
        for (size_t i = 0; i < etaX_.size(); i++) {
            cout << "abstraction " << i << ": " << *tau_[i] << '\n';
        }
    }


};


}

#endif /* ADAPTIVE_HH_ */
