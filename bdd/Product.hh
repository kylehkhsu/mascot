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
 *  \brief A class (derived from base Adaptive) that synthesizes multi-layer abstractions by taking products of a base system transition relation and auxicate transition relations.
 */

class Product: public Adaptive {
public:
    Cudd* prodDdmgr_;
    System* base_;
    vector<System*> auxs_;

    vector<double*> baseEtaXs_;
    vector<vector<double*>*> auxsEtaXs_;
    vector<double*> allTau_;

    vector<OdeSolver*> baseSolvers_;
    vector<vector<OdeSolver*>*> auxsSolvers_;

    vector<SymbolicSet*> baseXs_;
    vector<vector<SymbolicSet*>*> auxsXs_;
    vector<SymbolicSet*> baseX2s_;
    vector<vector<SymbolicSet*>*> auxsX2s_;

    SymbolicSet* baseU_;
    SymbolicSet* auxU_;

    vector<SymbolicSet*> baseTs_;
    vector<vector<SymbolicSet*>*> auxsTs_;

    vector<SymbolicSet*> prodTs_;

    System* system_;

    Product(char* logFile)
        : Adaptive(logFile) {}

    ~Product() {
        deleteVecArray(baseEtaXs_);
        deleteVecVecArray(auxsEtaXs_);
        deleteVec(allTau_);
        deleteVec(baseSolvers_);
        deleteVecVec(auxsSolvers_);
        deleteVec(baseXs_);
        deleteVecVec(auxsXs_);
        deleteVec(baseX2s_);
        deleteVecVec(auxsX2s_);
        delete baseU_;
        delete auxU_;
        deleteVec(baseTs_);
        deleteVecVec(auxsTs_);
        deleteVec(prodTs_);
        delete system_;
        delete prodDdmgr_;
    }

    void computeProducts() {
        vector<SymbolicSet*> curTs;


        for (int i = 0; i < *base_->numAbs_; i++) {
            SymbolicSet* T = new SymbolicSet(*baseTs_[i]);
            T->symbolicSet_ = baseTs_[i]->symbolicSet_;
            curTs.push_back(T);
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* curT = new SymbolicSet(*(curTs.back()), *((*auxsTs_[iAux])[i]));
                curT->symbolicSet_ = (curTs.back())->symbolicSet_ & ((*auxsTs_[iAux])[i])->symbolicSet_;
                curTs.push_back(curT);
                curT->printInfo(1);
            }

            SymbolicSet* prodT = new SymbolicSet(*(curTs.back()));
            prodT->symbolicSet_ = (curTs.back())->symbolicSet_;
            prodTs_.push_back(prodT);
        }
        deleteVec(curTs);

        checkMakeDir("T");
        saveVec(prodTs_, "T/T");
        clog << "Wrote prodTs_ to file.\n";
    }


    void initializeProduct(System* base, vector<System*> auxs) {
        prodDdmgr_ = new Cudd;
        base_ = base;
        auxs_ = auxs;

        if (auxs_.size() > 1) {
            error("Product currently supports only one auxicate.\n");
        }

        initializeProductEtaTau();
        initializeProductSolvers();
        initializeProductBDDs();

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            vector<SymbolicSet*>* auxTs = new vector<SymbolicSet*>;
            auxsTs_.push_back(auxTs);
        }

//        baseXs_[0]->printInfo(1);
//        (*auxsXs_[0])[0]->printInfo(1);

    }

    template<class O_type>
    void initializeAdaptive(int readXX, int readAbs, O_type addO) {
        initializeProductSystem();
        this->initialize(system_, readXX, readAbs, addO);
    }

    template<class sys_type, class rad_type, class x_type, class u_type>
    void computeBaseAbstractions(sys_type sysNext, rad_type radNext, x_type x, u_type u) {
        for (int i = 0; i < *base_->numAbs_; i++) {
            SymbolicModelGrowthBound<x_type, u_type> baseAb(baseXs_[i], baseU_, baseX2s_[i]);
            baseAb.computeTransitionRelation(sysNext, radNext, *baseSolvers_[i]);

            SymbolicSet T = baseAb.getTransitionRelation();
            SymbolicSet* baseT = new SymbolicSet(T);
            baseT->symbolicSet_ = T.symbolicSet_;

            baseTs_.push_back(baseT);

//            T.printInfo(1);
        }
    }

    template<class sys_type, class rad_type, class x_type, class u_type>
    void computeAuxAbstractions(sys_type sysNext, rad_type radNext, x_type x, u_type u, int iAux) {
        for (int i = 0; i < *base_->numAbs_; i++) {
            SymbolicModelGrowthBound<x_type, u_type> auxAb((*auxsXs_[iAux])[i], auxU_, (*auxsX2s_[iAux])[i]);
            auxAb.computeTransitionRelation(sysNext, radNext, *(*auxsSolvers_[iAux])[i]);

            SymbolicSet T = auxAb.getTransitionRelation();
            SymbolicSet* auxT = new SymbolicSet(*(*auxsXs_[iAux])[i], *(*auxsX2s_[iAux])[i]);
            auxT->symbolicSet_ = T.symbolicSet_.ExistAbstract(auxU_->getCube());

            auxsTs_[iAux]->push_back(auxT);

//            T.printInfo(1);
//            auxT->printInfo(1);
        }
    }

    void initializeProductSystem() {
        int dimX = 0;
        dimX = dimX + *base_->dimX_;
        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            System* aux = auxs_[iAux];
            dimX = dimX + *aux->dimX_;
        }
        double lbX[dimX] = {0};
        double ubX[dimX] = {0};
        double etaX[dimX] = {0};
        double etaRatio[dimX] = {0};

        for (int i = 0; i < *base_->dimX_; i++) {
            lbX[i] = base_->lbX_[i];
            ubX[i] = base_->ubX_[i];
            etaX[i] = base_->etaX_[i];
            etaRatio[i] = base_->etaRatio_[i];
        }

        int ind = *base_->dimX_;

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            System* aux = auxs_[iAux];
            for (int i = 0; i < *aux->dimX_; i++) {
                lbX[ind] = aux->lbX_[i];
                ubX[ind] = aux->ubX_[i];
                etaX[ind] = aux->etaX_[i];
                etaRatio[ind] = aux->etaRatio_[i];
                ind++;
            }
        }

        cout << "ind: " << ind << '\n';

        system_ = new System(dimX, lbX, ubX, etaX, *base_->tau_,
                                    *base_->dimU_, base_->lbU_, base_->ubU_, base_->etaU_,
                                    etaRatio, *base_->tauRatio_, *base_->nSubInt_, *base_->numAbs_);
    }

    void initializeProductSolvers() {
        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            vector<OdeSolver*>* auxSolvers = new vector<OdeSolver*>;
            auxsSolvers_.push_back(auxSolvers);
        }

        for (int i = 0; i < *base_->numAbs_; i++) {
            OdeSolver* baseSolver = new OdeSolver(*base_->dimX_, *base_->nSubInt_, allTau_[i][0]);
            baseSolvers_.push_back(baseSolver);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                System* aux = auxs_[iAux];
                OdeSolver* auxSolver = new OdeSolver(*aux->dimX_, *aux->nSubInt_, allTau_[i][0]);
                auxsSolvers_[iAux]->push_back(auxSolver);
            }
        }

        cout << "Initialized base's, auxs' ODE solvers.\n";
    }

    void initializeProductEtaTau() {

        double* baseEtaCur = new double[*base_->dimX_];
        for (int i = 0; i < *base_->dimX_; i++) {
            baseEtaCur[i] = base_->etaX_[i];
        }
        for (int i = 0; i < *base_->numAbs_; i++) {
            double* baseEtaX = new double[*base_->dimX_];
            for (int j = 0; j < *base_->dimX_; j++) {
                baseEtaX[j] = baseEtaCur[j];
            }
            baseEtaXs_.push_back(baseEtaX);
            for (int j = 0; j < *base_->dimX_; j++) {
                baseEtaCur[j] /= base_->etaRatio_[j];
            }
        }
        delete[] baseEtaCur;

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            System* aux = auxs_[iAux];
            vector<double*>* auxEtaX = new vector<double*>;

            double* auxEtaCur = new double[*aux->dimX_];
            for (int i = 0; i < *aux->dimX_; i++) {
                auxEtaCur[i] = aux->etaX_[i];
            }
            for (int i = 0; i < *aux->dimX_; i++) {
                double* auxEtai = new double[*aux->dimX_];
                for (int j = 0; j < *aux->dimX_; j++) {
                    auxEtai[j] = auxEtaCur[j];
                }
                auxEtaX->push_back(auxEtai);
                for (int j = 0; j < *aux->dimX_; j++) {
                    auxEtaCur[j] /= aux->etaRatio_[j];
                }
            }
            delete[] auxEtaCur;
            auxsEtaXs_.push_back(auxEtaX);
        }

        double* allTauCur = new double;
        *allTauCur = *base_->tau_;

        for (int i = 0; i < *base_->numAbs_; i++) {
            double* allTaui = new double;
            *allTaui = *allTauCur;
            allTau_.push_back(allTaui);
            *allTauCur /= *base_->tauRatio_;
        }
        delete allTauCur;

        clog << "Initialized base's, auxs' etaXs, taus.\n";

    }

    void initializeProductBDDs() {
        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            vector<SymbolicSet*>* auxXs = new vector<SymbolicSet*>;
            auxsXs_.push_back(auxXs);
            vector<SymbolicSet*>* auxX2s = new vector<SymbolicSet*>;
            auxsX2s_.push_back(auxX2s);
        }

        for (int i = 0; i < *base_->numAbs_; i++) {
            SymbolicSet* baseX = new SymbolicSet(*prodDdmgr_, *base_->dimX_, base_->lbX_, base_->ubX_, baseEtaXs_[i], allTau_[i][0]);
            baseX->addGridPoints();
            baseXs_.push_back(baseX);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                System* aux = auxs_[iAux];
                SymbolicSet* auxX = new SymbolicSet(*prodDdmgr_, *aux->dimX_, aux->lbX_, aux->ubX_, (*auxsEtaXs_[iAux])[i], allTau_[i][0]);
                auxX->addGridPoints();
                auxsXs_[iAux]->push_back(auxX);
            }
        }

        for (int i = 0; i < *base_->numAbs_; i++) {
            SymbolicSet* baseX2 = new SymbolicSet(*baseXs_[i], 1);
            baseX2s_.push_back(baseX2);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* auxX2 = new SymbolicSet(*((*auxsXs_[iAux])[i]), 1);
                auxsX2s_[iAux]->push_back(auxX2);
            }
        }

        baseU_ = new SymbolicSet(*prodDdmgr_, *base_->dimU_, base_->lbU_, base_->ubU_, base_->etaU_, 0);
        baseU_->addGridPoints();

        System* aux = auxs_[0];
        auxU_ = new SymbolicSet(*prodDdmgr_, *aux->dimU_, aux->lbU_, aux->ubU_, aux->etaU_, 0);
        auxU_->addGridPoints();

        clog << "Initialized base's, auxs' Xs, X2s, U.\n";
    }

};


}

#endif /* Product_HH_ */
