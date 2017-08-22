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
    vector<SymbolicSet*> baseX2s_;
    vector<SymbolicSet*> baseXXs_;

    vector<vector<SymbolicSet*>*> auxsXs_;
    vector<vector<SymbolicSet*>*> auxsX2s_;
    vector<vector<SymbolicSet*>*> auxsXXs_;

    vector<vector<SymbolicSet*>*> prodsXs_;
    vector<vector<SymbolicSet*>*> prodsX2s_;
    vector<vector<SymbolicSet*>*> prodsXXs_;

    vector<vector<SymbolicSet*>*> prodsCs_;
    vector<vector<SymbolicSet*>*> prodsTs_;

    SymbolicSet* baseU_;
    SymbolicSet* auxU_;

    vector<SymbolicSet*> baseTs_;
    vector<vector<SymbolicSet*>*> auxsTs_;

    vector<SymbolicSet*> prodTs_;

    vector<BDD*> baseCubesX_;
    vector<BDD*> baseCubesX2_;
    vector<vector<BDD*>*> auxsCubesX_;
    vector<vector<BDD*>*> auxsCubesX2_;
    vector<vector<BDD*>*> prodsCubesX_; // useless?
    vector<vector<BDD*>*> prodsCubesX2_; // useless?
    vector<BDD*> baseNotXUVars_;
    vector<BDD*> baseNotXVars_;
    vector<vector<BDD*>*> prodsNotXUVars_;
    vector<vector<BDD*>*> prodsNotXVars_;

    vector<vector<int*>*> prodsPermutesXtoX2_;
    vector<vector<int*>*> prodsPermutesX2toX_;

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
        deleteVec(baseX2s_);
        deleteVec(baseXXs_);
        deleteVecVec(auxsXs_);
        deleteVecVec(auxsX2s_);
        deleteVecVec(auxsXXs_);
        deleteVecVec(prodsXs_);
        deleteVecVec(prodsX2s_);
        deleteVecVec(prodsXXs_);

        delete baseU_;
        delete auxU_;

        deleteVecVec(prodsCs_);
        deleteVecVec(prodsTs_);

        deleteVec(baseTs_);
        deleteVecVec(auxsTs_);
        deleteVec(prodTs_);

        deleteVec(baseCubesX_);
        deleteVec(baseCubesX2_);
        deleteVecVec(auxsCubesX_);
        deleteVecVec(auxsCubesX2_);
        deleteVecVec(prodsCubesX_);
        deleteVecVec(prodsCubesX2_);
        deleteVec(baseNotXUVars_);
        deleteVec(baseNotXVars_);
        deleteVecVec(prodsNotXUVars_);
        deleteVecVec(prodsNotXVars_);

        deleteVecVecArray(prodsPermutesXtoX2_);
        deleteVecVecArray(prodsPermutesX2toX_);

        delete system_;
        delete prodDdmgr_;
    }

    void constructAbstractions() {
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                ((*prodsTs_[iAux])[iAbs])->symbolicSet_ = baseTs_[iAbs]->symbolicSet_ & ((*auxsTs_[iAux])[iAbs])->symbolicSet_;
            }
        }

        printVecVec(prodsTs_, "prodsTs");
    }

    void computeProducts() {
        vector<SymbolicSet*> curTs;


        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicSet* T = new SymbolicSet(*baseTs_[iAbs]);
            T->symbolicSet_ = baseTs_[iAbs]->symbolicSet_;
            curTs.push_back(T);
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* curT = new SymbolicSet(*(curTs.back()), *((*auxsTs_[iAux])[iAbs]));
                curT->symbolicSet_ = (curTs.back())->symbolicSet_ & ((*auxsTs_[iAux])[iAbs])->symbolicSet_;
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
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicModelGrowthBound<x_type, u_type> baseAb(baseXs_[iAbs], baseU_, baseX2s_[iAbs]);
            baseAb.computeTransitionRelation(sysNext, radNext, *baseSolvers_[iAbs]);

            SymbolicSet T = baseAb.getTransitionRelation();
            SymbolicSet* baseT = new SymbolicSet(T);
            baseT->symbolicSet_ = T.symbolicSet_;

            baseTs_.push_back(baseT);

//            T.printInfo(1);
        }
    }

    template<class sys_type, class rad_type, class x_type, class u_type>
    void computeAuxAbstractions(sys_type sysNext, rad_type radNext, x_type x, u_type u, size_t iAux) {
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicModelGrowthBound<x_type, u_type> auxAb((*auxsXs_[iAux])[iAbs], auxU_, (*auxsX2s_[iAux])[iAbs]);
            auxAb.computeTransitionRelation(sysNext, radNext, *(*auxsSolvers_[iAux])[iAbs]);

            SymbolicSet T = auxAb.getTransitionRelation();
            SymbolicSet* auxT = new SymbolicSet(*(*auxsXs_[iAux])[iAbs], *(*auxsX2s_[iAux])[iAbs]);
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

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            OdeSolver* baseSolver = new OdeSolver(*base_->dimX_, *base_->nSubInt_, allTau_[iAbs][0]);
            baseSolvers_.push_back(baseSolver);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                System* aux = auxs_[iAux];
                OdeSolver* auxSolver = new OdeSolver(*aux->dimX_, *base_->nSubInt_, allTau_[iAbs][0]);
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
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
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
            vector<SymbolicSet*>* auxXXs = new vector<SymbolicSet*>;
            auxsXXs_.push_back(auxXXs);
            vector<SymbolicSet*>* prodXs = new vector<SymbolicSet*>;
            prodsXs_.push_back(prodXs);
            vector<SymbolicSet*>* prodX2s = new vector<SymbolicSet*>;
            prodsX2s_.push_back(prodX2s);
            vector<SymbolicSet*>* prodXXs = new vector<SymbolicSet*>;
            prodsXXs_.push_back(prodXXs);
            vector<SymbolicSet*>* prodCs = new vector<SymbolicSet*>;
            prodsCs_.push_back(prodCs);
            vector<SymbolicSet*>* prodTs = new vector<SymbolicSet*>;
            prodsTs_.push_back(prodTs);

            vector<BDD*>* auxCubesX = new vector<BDD*>;
            vector<BDD*>* auxCubesX2 = new vector<BDD*>;
            vector<BDD*>* prodCubesX = new vector<BDD*>;
            vector<BDD*>* prodCubesX2 = new vector<BDD*>;
            auxsCubesX_.push_back(auxCubesX);
            auxsCubesX2_.push_back(auxCubesX2);
            prodsCubesX_.push_back(prodCubesX);
            prodsCubesX2_.push_back(prodCubesX2);

            vector<int*>* prodPermutesXtoX2 = new vector<int*>;
            vector<int*>* prodPermutesX2toX = new vector<int*>;
            prodsPermutesXtoX2_.push_back(prodPermutesXtoX2);
            prodsPermutesX2toX_.push_back(prodPermutesX2toX);

            vector<BDD*>* prodNotXUVars = new vector<BDD*>;
            vector<BDD*>* prodNotXVars = new vector<BDD*>;
            prodsNotXUVars_.push_back(prodNotXUVars);
            prodsNotXVars_.push_back(prodNotXVars);
        }

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicSet* baseX = new SymbolicSet(*prodDdmgr_, *base_->dimX_, base_->lbX_, base_->ubX_, baseEtaXs_[iAbs], allTau_[iAbs][0]);
            baseX->addGridPoints();
            baseXs_.push_back(baseX);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                System* aux = auxs_[iAux];
                SymbolicSet* auxX = new SymbolicSet(*prodDdmgr_, *aux->dimX_, aux->lbX_, aux->ubX_, (*auxsEtaXs_[iAux])[iAbs], allTau_[iAbs][0]);
                auxX->addGridPoints();
                auxsXs_[iAux]->push_back(auxX);

                SymbolicSet* prodX = new SymbolicSet(*baseXs_[iAbs], *auxX);
                prodX->addGridPoints();
                prodsXs_[iAux]->push_back(prodX);
            }
        }

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicSet* baseX2 = new SymbolicSet(*baseXs_[iAbs], 1);
            baseX2s_.push_back(baseX2);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* auxX2 = new SymbolicSet(*((*auxsXs_[iAux])[iAbs]), 1);
                auxsX2s_[iAux]->push_back(auxX2);

                SymbolicSet* prodX2 = new SymbolicSet(*baseX2s_[iAbs], *((*auxsX2s_[iAux])[iAbs]));
                prodsX2s_[iAux]->push_back(prodX2);
            }
        }
        clog << "Initialized base's, auxs', prods' Xs, X2s to full.\n";

        for (int iAbs = 0; iAbs < *base_->numAbs_ - 1; iAbs++) {
            SymbolicSet* baseXX = new SymbolicSet(*baseXs_[iAbs], *baseXs_[iAbs+1]);
            this->mapAbstractions(baseXs_[iAbs], baseXs_[iAbs+1], baseXX, iAbs);
            baseXXs_.push_back(baseXX);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* auxXX = new SymbolicSet(*((*auxsXs_[iAux])[iAbs]), *((*auxsXs_[iAux])[iAbs+1]));
                this->mapAbstractions((*auxsXs_[iAux])[iAbs], (*auxsXs_[iAux])[iAbs+1], auxXX, iAbs);
                auxsXXs_[iAux]->push_back(auxXX);

                SymbolicSet* prodXX = new SymbolicSet(*((*prodsXs_[iAux])[iAbs]), *((*prodsXs_[iAux])[iAbs+1]));
                prodXX->symbolicSet_ = baseXX->symbolicSet_ & auxXX->symbolicSet_;
                prodsXXs_[iAux]->push_back(prodXX);
            }
        }
        clog << "Initialized base's, auxs', prods' XXs via mapAbstractions.\n";

        baseU_ = new SymbolicSet(*prodDdmgr_, *base_->dimU_, base_->lbU_, base_->ubU_, base_->etaU_, 0);
        baseU_->addGridPoints();

        System* aux = auxs_[0];
        auxU_ = new SymbolicSet(*prodDdmgr_, *aux->dimU_, aux->lbU_, aux->ubU_, aux->etaU_, 0);
        auxU_->addGridPoints();

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* prodC = new SymbolicSet(*((*prodsXs_[iAux])[iAbs]), *baseU_);
                prodsCs_[iAux]->push_back(prodC);
            }
        }
        clog << "Initialized prods' Cs to empty.\n";

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* prodT = new SymbolicSet(*((*prodsCs_[iAux])[iAbs]), *((*prodsX2s_[iAux])[iAbs]));
                prodsTs_[iAux]->push_back(prodT);
            }
        }
        clog << "Initialized prods' Ts to empty.\n";

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            BDD* baseCubeX = new BDD;
            BDD* baseCubeX2 = new BDD;
            *baseCubeX = baseXs_[iAbs]->getCube();
            *baseCubeX2 = baseX2s_[iAbs]->getCube();
            baseCubesX_.push_back(baseCubeX);
            baseCubesX2_.push_back(baseCubeX2);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                BDD* auxCubeX = new BDD;
                BDD* auxCubeX2 = new BDD;
                BDD* prodCubeX = new BDD;
                BDD* prodCubeX2 = new BDD;

                *auxCubeX = ((*auxsXs_[iAux])[iAbs])->getCube();
                *auxCubeX2 = ((*auxsX2s_[iAux])[iAbs])->getCube();
                *prodCubeX = ((*prodsXs_[iAux])[iAbs])->getCube();
                *prodCubeX2 = ((*prodsX2s_[iAux])[iAbs])->getCube();

                auxsCubesX_[iAux]->push_back(auxCubeX);
                auxsCubesX2_[iAux]->push_back(auxCubeX2);
                prodsCubesX_[iAux]->push_back(prodCubeX);
                prodsCubesX2_[iAux]->push_back(prodCubeX2);
            }
        }
        clog << "Initialized base's, auxs', prods' cubesX, cubesX2.\n";

        int numBDDVars = 0;
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            numBDDVars += baseXs_[iAbs]->nvars_;
            numBDDVars += baseX2s_[iAbs]->nvars_;
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                numBDDVars += ((*auxsXs_[iAux])[iAbs])->nvars_;
                numBDDVars += ((*auxsX2s_[iAux])[iAbs])->nvars_;
            }
        }
        numBDDVars += baseU_->nvars_;
        numBDDVars += auxU_->nvars_;
        clog << "Number of BDD variables: " << numBDDVars << '\n';

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                int* prodPermuteXtoX2 = new int[numBDDVars];
                int* prodPermuteX2toX = new int[numBDDVars];
                for (int i = 0; i < numBDDVars; i++) {
                    prodPermuteXtoX2[i] = 0;
                    prodPermuteX2toX[i] = 0;
                }
                for (size_t iVar = 0; iVar < ((*prodsXs_[iAux])[iAbs])->nvars_; iVar++) {
                    prodPermuteXtoX2[((*prodsXs_[iAux])[iAbs])->idBddVars_[iVar]] = ((*prodsX2s_[iAux])[iAbs])->idBddVars_[iVar];
                    prodPermuteX2toX[((*prodsX2s_[iAux])[iAbs])->idBddVars_[iVar]] = ((*prodsXs_[iAux])[iAbs])->idBddVars_[iVar];
                }
                prodsPermutesXtoX2_[iAux]->push_back(prodPermuteXtoX2);
                prodsPermutesX2toX_[iAux]->push_back(prodPermuteX2toX);
            }
        }
        clog << "Initialized prods' permutes.\n";

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                BDD* prodNotXUVar = new BDD;
                BDD* prodNotXVar = new BDD;
                *prodNotXUVar = prodDdmgr_->bddOne();

                for (int jAbs = 0; jAbs < *base_->numAbs_; jAbs++) {
                    *prodNotXUVar &= *baseCubesX2_[jAbs];
                    if (!(iAbs == jAbs)) {
                        *prodNotXUVar &= *baseCubesX_[jAbs];
                    }
                    for (size_t jAux = 0; jAux < auxs_.size(); jAux++) {
                        *prodNotXUVar &= *((*auxsCubesX2_[jAux])[jAbs]);
                        if (!((iAbs == jAbs) && (iAux == jAux))) {
                            *prodNotXUVar &= *((*auxsCubesX_[jAux])[jAbs]);
                        }
                    }
                }

                *prodNotXUVar &= auxU_->getCube();
                *prodNotXVar = *prodNotXUVar & baseU_->getCube();
                prodsNotXUVars_[iAux]->push_back(prodNotXUVar);
                prodsNotXVars_[iAux]->push_back(prodNotXVar);
            }

            BDD* baseNotXUVar = new BDD;
            BDD* baseNotXVar = new BDD;
            *baseNotXUVar = prodDdmgr_->bddOne();

            *baseNotXUVar &= *((*prodsNotXUVars_[0])[iAbs]);
            *baseNotXUVar &= *((*auxsCubesX_[0])[iAbs]);
            *baseNotXVar = *baseNotXUVar & baseU_->getCube();
            baseNotXUVars_.push_back(baseNotXUVar);
            baseNotXVars_.push_back(baseNotXVar);
        }

        printVec(baseXs_, "baseXs");
//        printVec(baseXXs_, "baseXXs");

//        printVecVec(auxsXs_, "auxsXs");
//        printVecVec(auxsXXs_, "auxsXXs");

        printVecVec(prodsXs_, "prodsXs");
//        printVecVec(prodsXXs_, "prodsXXs");

//        printVecVec(prodsTs_, "prodsTs");


        checkNotVars();




    }

    void checkNotVars() {
        // for 2A1P
        SymbolicSet d1(*((*prodsXs_[0])[0]), *((*prodsXs_[0])[1]));
        SymbolicSet d2(d1, *((*prodsX2s_[0])[0]));
        SymbolicSet d3(d2, *((*prodsX2s_[0])[1]));
        SymbolicSet d4(d3, *baseU_);
        SymbolicSet all(d4, *auxU_);
//        all.symbolicSet_ = *((*auxsCubesX_[0])[0]);
//        all.printInfo(2);

//        all.symbolicSet_ = *((*prodsNotXVars_[0])[0]);
//        all.printInfo(2);

        all.symbolicSet_ = *baseNotXVars_[1];
        cout << "baseNotXVars_[1]:\n";
        all.printInfo(2);
    }

};


}

#endif /* Product_HH_ */
