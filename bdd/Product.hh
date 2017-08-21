/*! \file Product.hh
 *  Contains the Product class. */

#ifndef PRODUCT_HH_
#define PRODUCT_HH_

#include <vector>

#include "Adaptive.hh"
#include "System.hh"
#include "SymbolicModelGrowthBound.hh"
#include "RungeKutta4.hh"

using namespace helper;
using std::vector;
using std::cout;

namespace scots {

/*! \class Product
 *  \brief A class (derived from base Adaptive) that synthesizes multi-layer abstractions by taking products of a base system transition relation and predicate transition relations.
 */
template<class X_type, class U_type>
class Product: public Adaptive<X_type, U_type> {
public:
    Cudd* ddmgr_;
    System* dyn_;
    vector<System*> preds_;

    vector<double*> dynEtaXs_;
    vector<vector<double*>*> predsEtaXs_;
    vector<double*> allTau_;

    vector<OdeSolver*> dynSolvers_;
    vector<vector<OdeSolver*>*> predsSolvers_;

    vector<SymbolicSet*> dynXs_;
    vector<vector<SymbolicSet*>*> predsXs_;
    vector<SymbolicSet*> dynX2s_;
    vector<vector<SymbolicSet*>*> predsX2s_;

    SymbolicSet* dynU_;
    SymbolicSet* predU_;

    vector<SymbolicSet*> dynTs_;
    vector<vector<SymbolicSet*>*> predsTs_;

    Product(char* logFile)
        : Adaptive<X_type, U_type>(logFile) {}

    ~Product() {
        deleteVecArray(dynEtaXs_);
        deleteVecVecArray(predsEtaXs_);
        deleteVec(allTau_);
        deleteVec(dynSolvers_);
        deleteVecVec(predsSolvers_);
        deleteVec(dynXs_);
        deleteVecVec(predsXs_);
        deleteVec(dynX2s_);
        deleteVecVec(predsX2s_);
        delete dynU_;
        delete predU_;
        deleteVec(dynTs_);
        deleteVecVec(predsTs_);
        delete ddmgr_;
    }

    void initializeProduct(System* dyn, vector<System*> preds) {
        ddmgr_ = new Cudd;
        dyn_ = dyn;
        preds_ = preds;

        if (preds_.size() > 1) {
            error("Product currently supports only one predicate.\n");
        }

        initializeProductEtaTau();
        initializeProductSolvers();
        initializeProductBDDs();

        for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
            vector<SymbolicSet*>* predTs = new vector<SymbolicSet*>;
            predsTs_.push_back(predTs);
        }

//        dynXs_[0]->printInfo(1);
//        (*predsXs_[0])[0]->printInfo(1);

    }




    template<class sys_type, class rad_type, class x_type, class u_type>
    void computeDynAbstractions(sys_type sysNext, rad_type radNext, x_type x, u_type u) {
        for (int i = 0; i < *dyn_->numAbs_; i++) {
            SymbolicModelGrowthBound<x_type, u_type> dynAb(dynXs_[i], dynU_, dynX2s_[i]);
            dynAb.computeTransitionRelation(sysNext, radNext, *dynSolvers_[i]);

            SymbolicSet T = dynAb.getTransitionRelation();
            SymbolicSet* dynT = new SymbolicSet(T);
            dynT->symbolicSet_ = T.symbolicSet_;

            dynTs_.push_back(dynT);

            T.printInfo(1);
        }
    }

    template<class sys_type, class rad_type, class x_type, class u_type>
    void computePredAbstractions(sys_type sysNext, rad_type radNext, x_type x, u_type u, int iPred) {
        for (int i = 0; i < *dyn_->numAbs_; i++) {
            SymbolicModelGrowthBound<x_type, u_type> predAb((*predsXs_[iPred])[i], predU_, (*predsX2s_[iPred])[i]);
            predAb.computeTransitionRelation(sysNext, radNext, *(*predsSolvers_[iPred])[i]);

            SymbolicSet T = predAb.getTransitionRelation();
            SymbolicSet* predT = new SymbolicSet(T);
            predT->symbolicSet_ = T.symbolicSet_;

            predsTs_[iPred]->push_back(predT);

            T.printInfo(1);
        }
    }

    System* getProductSystem() {
        int dimX = 0;
        dimX = dimX + *dyn_->dimX_;
        for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
            System* pred = preds_[iPred];
            dimX = dimX + *pred->dimX_;
        }
        double lbX[dimX] = {0};
        double ubX[dimX] = {0};
        double etaX[dimX] = {0};
        double etaRatio[dimX] = {0};

        for (int i = 0; i < *dyn_->dimX_; i++) {
            lbX[i] = dyn_->lbX_[i];
            ubX[i] = dyn_->ubX_[i];
            etaX[i] = dyn_->etaX_[i];
            etaRatio[i] = dyn_->etaRatio_[i];
        }

        int ind = *dyn_->dimX_;

        for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
            System* pred = preds_[iPred];
            for (int i = 0; i < *pred->dimX_; i++) {
                lbX[ind] = pred->lbX_[i];
                ubX[ind] = pred->ubX_[i];
                etaX[ind] = pred->etaX_[i];
                etaRatio[ind] = pred->etaRatio_[i];
                ind++;
            }
        }

        cout << "ind: " << ind << '\n';

        System* system = new System(dimX, lbX, ubX, etaX, *dyn_->tau_,
                                    *dyn_->dimU_, dyn_->lbU_, dyn_->ubU_, dyn_->etaU_,
                                    etaRatio, *dyn_->tauRatio_, *dyn_->nSubInt_, *dyn_->numAbs_);
        return system;
    }

    void initializeProductSolvers() {
        for (int i = 0; i < *dyn_->numAbs_; i++) {
            OdeSolver* dynSolver = new OdeSolver(*dyn_->dimX_, *dyn_->nSubInt_, allTau_[i][0]);
            dynSolvers_.push_back(dynSolver);

            vector<OdeSolver*>* predSolvers = new vector<OdeSolver*>;
            for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
                System* pred = preds_[iPred];
                OdeSolver* predSolver = new OdeSolver(*pred->dimX_, *pred->nSubInt_, allTau_[i][0]);
                predSolvers->push_back(predSolver);
            }
            predsSolvers_.push_back(predSolvers);
        }

        cout << "Initialized dyn's, preds' ODE solvers.\n";
    }

    void initializeProductEtaTau() {

        double* dynEtaCur = new double[*dyn_->dimX_];
        for (int i = 0; i < *dyn_->dimX_; i++) {
            dynEtaCur[i] = dyn_->etaX_[i];
        }
        for (int i = 0; i < *dyn_->numAbs_; i++) {
            double* dynEtaX = new double[*dyn_->dimX_];
            for (int j = 0; j < *dyn_->dimX_; j++) {
                dynEtaX[j] = dynEtaCur[j];
            }
            dynEtaXs_.push_back(dynEtaX);
            for (int j = 0; j < *dyn_->dimX_; j++) {
                dynEtaCur[j] /= dyn_->etaRatio_[j];
            }
        }
        delete[] dynEtaCur;

        for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
            System* pred = preds_[iPred];
            vector<double*>* predEtaX = new vector<double*>;

            double* predEtaCur = new double[*pred->dimX_];
            for (int i = 0; i < *pred->dimX_; i++) {
                predEtaCur[i] = pred->etaX_[i];
            }
            for (int i = 0; i < *pred->dimX_; i++) {
                double* predEtai = new double[*pred->dimX_];
                for (int j = 0; j < *pred->dimX_; j++) {
                    predEtai[j] = predEtaCur[j];
                }
                predEtaX->push_back(predEtai);
                for (int j = 0; j < *pred->dimX_; j++) {
                    predEtaCur[j] /= pred->etaRatio_[j];
                }
            }
            delete[] predEtaCur;
            predsEtaXs_.push_back(predEtaX);
        }

        double* allTauCur = new double;
        *allTauCur = *dyn_->tau_;

        for (int i = 0; i < *dyn_->numAbs_; i++) {
            double* allTaui = new double;
            *allTaui = *allTauCur;
            allTau_.push_back(allTaui);
            *allTauCur /= *dyn_->tauRatio_;
        }
        delete allTauCur;

        clog << "Initialized dyn's, preds' etaXs, taus.\n";

    }

    void initializeProductBDDs() {
        for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
            vector<SymbolicSet*>* predXs = new vector<SymbolicSet*>;
            predsXs_.push_back(predXs);
            vector<SymbolicSet*>* predX2s = new vector<SymbolicSet*>;
            predsX2s_.push_back(predX2s);
        }

        for (int i = 0; i < *dyn_->numAbs_; i++) {
            SymbolicSet* dynX = new SymbolicSet(*ddmgr_, *dyn_->dimX_, dyn_->lbX_, dyn_->ubX_, dynEtaXs_[i], allTau_[i][0]);
            dynX->addGridPoints();
            dynXs_.push_back(dynX);

            for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
                System* pred = preds_[iPred];
                SymbolicSet* predX = new SymbolicSet(*ddmgr_, *pred->dimX_, pred->lbX_, pred->ubX_, (*predsEtaXs_[iPred])[i], allTau_[i][0]);
                predX->addGridPoints();
                predsXs_[iPred]->push_back(predX);
            }
        }

        for (int i = 0; i < *dyn_->numAbs_; i++) {
            SymbolicSet* dynX2 = new SymbolicSet(*dynXs_[i], 1);
            dynX2s_.push_back(dynX2);

            for (size_t iPred = 0; iPred < preds_.size(); iPred++) {
                SymbolicSet* predX2 = new SymbolicSet(*((*predsXs_[iPred])[i]), 1);
                predsX2s_[iPred]->push_back(predX2);
            }
        }

        dynU_ = new SymbolicSet(*ddmgr_, *dyn_->dimU_, dyn_->lbU_, dyn_->ubU_, dyn_->etaU_, 0);
        dynU_->addGridPoints();

        System* pred = preds_[0];
        predU_ = new SymbolicSet(*ddmgr_, *pred->dimU_, pred->lbU_, pred->ubU_, pred->etaU_, 0);
        predU_->addGridPoints();

        clog << "Initialized dyn's, preds' Xs, X2s, U.\n";
    }

};




}

#endif /* Product_HH_ */
