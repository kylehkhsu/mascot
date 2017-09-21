/*! \file Composition.hh
 *  Contains the Composition class. */

#ifndef Composition_HH_
#define Composition_HH_

#include <vector>

#include "Helper.hh"
#include "System.hh"
#include "SymbolicModelGrowthBound.hh"
#include "RungeKutta4.hh"

using std::vector;
using std::cout;

namespace scots {
/*! \class Composition
 *  \brief A class that constructs multi-layer abstractions of compositional systems.
 */
class Composition {
public:
    Cudd* ddmgr_; /*!< Decision diagram manager. */
    System* base_; /*!< Base system. */
    vector<System*> auxs_; /*!< Auxiliary systems.*/

    vector<double*> baseEtaXs_; /*!< Grid parameters for base system abstractions. */
    vector<vector<double*>*> auxsEtaXs_; /*!< Grid parameters for auxiliary system abstractions. */
    vector<vector<double*>*> prodsEtaXs_;
    vector<double*> prodsEtaRatio_;
    vector<double*> taus_; /*!< Time sampling parameters. */
    vector<OdeSolver*> baseSolvers_; /*!< Runge-Kutta solvers for base system abstractions. */
    vector<vector<OdeSolver*>*> auxsSolvers_; /*!< Runge-Kutta solvers for auxiliary system abstractions. */

    vector<SymbolicSet*> baseXs_; /*!< State spaces of base system abstractions. */
    vector<SymbolicSet*> baseX2s_; /*!< Post-state spaces of base system abstractions. */
    vector<SymbolicSet*> baseXXs_; /*!< Projection mappings between consecutive base system abstractions. */
    vector<vector<SymbolicSet*>*> auxsXs_; /*!< State spaces of auxiliary system abstractions. */
    vector<vector<SymbolicSet*>*> auxsX2s_; /*!< Post-state spaces of auxiliary system abstractions. */
    vector<vector<SymbolicSet*>*> auxsXXs_; /*!< Projection mappings between consecutive auxiliary system abstractions. */

    vector<vector<SymbolicSet*>*> prodsXs_; /*!< State spaces of product abstractions. */
    vector<vector<SymbolicSet*>*> prodsX2s_; /*!< Post-state spaces of product abstractions. */
    vector<vector<SymbolicSet*>*> prodsXXs_; /*!< Projection mappings between consecutive product abstractions. */

    vector<SymbolicSet*> baseOs_; /*!< Instance of baseXs_ containing obstacle states. */
    vector<vector<SymbolicSet*>*> prodsGs_; /*!< Instance of prodsXs_ containing goal states. */

    vector<vector<SymbolicSet*>*> prodsZs_; /*!< Instance of prodsXs_ containing controller domain states. */
    vector<vector<SymbolicSet*>*> prodsValidZs_; /*!< Instance of prodsXs_ containing controller domain states declared valid. */
    vector<vector<SymbolicSet*>*> prodsFinalZs_; /*!< Contains saved controller domain states. */

    vector<vector<SymbolicSet*>*> prodsYs_; /*!< Instance of prodsXs_ containing converged unified controller domain states. */
    vector<vector<SymbolicSet*>*> prodsPreYs_; /*!< Instance of prodsXs_ containing controllable predecessors of converged unified controller domain states. */
    vector<vector<SymbolicSet*>*> prodsPrevYs_; /*!< Instance of prodsXs_ containing previous iteration's prodsYs_. */
    vector<vector<SymbolicSet*>*> prodsPrevPreYs_; /*!< Instance of prodsXs_ containing previous iterations' prodsPreYs_. */

    vector<vector<SymbolicSet*>*> prodsCs_; /*!< Controllers of product abstractions. */
    vector<vector<SymbolicSet*>*> prodsValidCs_; /*!< Controllers of product declared valid. */
    vector<vector<SymbolicSet*>*> prodsFinalCs_; /*!< Contains saved controllers. Each vector is the sequence of controllers that brings the state to a particular target. Written to child directory 'C'. */

    vector<SymbolicSet*> baseTs_; /*!< Transition systems of base system. */
    vector<vector<SymbolicSet*>*> auxsTs_; /*!< Transition systems of auxiliary systems. */
    vector<vector<SymbolicSet*>*> prodsTs_; /*!< Transition systems of product systems. */
    vector<vector<SymbolicSet*>*> prodsTTs_; /*!< Transition systems of product systems with post-states existentially abstracted */

    SymbolicSet* baseU_; /*!< Input space of base system. */
    SymbolicSet* auxU_; /*!< Input space (dummy) of auxiliary systems. */

    vector<BDD*> baseCubesX_; /*!< Cubes for state spaces of base system abstractions. */
    vector<BDD*> baseCubesX2_; /*!< Cubes for post-state spaces of base system abstractions. */
    vector<vector<BDD*>*> auxsCubesX_; /*!< Cubes for state spaces of auxiliary systems' abstractions. */
    vector<vector<BDD*>*> auxsCubesX2_; /*!< Cubes for post-state spaces of auxiliary systems' abstractions. */
    vector<BDD*> baseNotXUVars_; /*!< Used for existential quantification for which result's domain is base system controller space. */
    vector<BDD*> baseNotXVars_; /*!< Used for existential quantification for which result's domain is base system state space. */
    vector<vector<BDD*>*> prodsNotXUVars_; /*!< Used for existential quantification for which result's domain is product system controller space. */
    vector<vector<BDD*>*> prodsNotXVars_; /*!< Used for existential quantification for which result's domain is product system state space. */
    vector<vector<int*>*> prodsPermutesXtoX2_; /*!< Used for projecting from state space to post-state space of product system abstractions. */
    vector<vector<int*>*> prodsPermutesX2toX_; /*!< Used for projecting from post-state space to state space of product system abstractions. */

    vector<vector<int*>*> prodsPermutesCoarser_;
    vector<vector<int*>*> prodsPermutesFiner_;
    vector<vector<BDD*>*> prodsCubesCoarser_;

    int* baseNumFiner_; /*!< Number of finer cells that compose a coarser cell for consecutive base system abstractions. */
    int* prodsNumFiner_; /*!< Number of finer cells that compose a coarser cell for consecutive product system abstractions. */

    /*! Constructor for a Composition object.
     *  \param[in]  logFile     Filename of program log.
     */
    Composition(char* logFile) {
        freopen(logFile, "w", stderr);
        clog << logFile << '\n';
    }

    /*! Destructor for a Composition object. */
    ~Composition() {
        deleteVecArray(baseEtaXs_);
        deleteVecVecArray(auxsEtaXs_);
        deleteVecVecArray(prodsEtaXs_);
        deleteVecArray(prodsEtaRatio_);
        deleteVec(taus_);
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

        deleteVec(baseOs_);
        deleteVecVec(prodsGs_);

        deleteVecVec(prodsZs_);
        deleteVecVec(prodsValidZs_);
        deleteVecVec(prodsFinalZs_);

        deleteVecVec(prodsYs_);
        deleteVecVec(prodsPreYs_);
        deleteVecVec(prodsPrevYs_);
        deleteVecVec(prodsPrevPreYs_);

        deleteVecVec(prodsCs_);
        deleteVecVec(prodsValidCs_);
        deleteVecVec(prodsFinalCs_);

        deleteVec(baseTs_);
        deleteVecVec(auxsTs_);
        deleteVecVec(prodsTs_);
        deleteVecVec(prodsTTs_);

        delete baseU_;
        delete auxU_;

        deleteVec(baseCubesX_);
        deleteVec(baseCubesX2_);
        deleteVecVec(auxsCubesX_);
        deleteVecVec(auxsCubesX2_);
        deleteVec(baseNotXUVars_);
        deleteVec(baseNotXVars_);
        deleteVecVec(prodsNotXUVars_);
        deleteVecVec(prodsNotXVars_);

        deleteVecVecArray(prodsPermutesXtoX2_);
        deleteVecVecArray(prodsPermutesX2toX_);

        deleteVecVecArray(prodsPermutesCoarser_);
        deleteVecVecArray(prodsPermutesFiner_);
        deleteVecVec(prodsCubesCoarser_);

        fclose(stderr);

        delete ddmgr_;
    }

    /*! Initializes goal set specification for a single product system for all layers in prodsG_[iAux].
     *  \param[in] addG     Pointer to function that adds goal states to a SymbolicSet.
     *  \param[in] iAux     0-index of product system.
     */
    template<class G_t>
    void initializeProdGs(G_t addG, int iAux) {
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicSet* prodG = new SymbolicSet(*((*prodsXs_[iAux])[iAbs]));
            addG(prodG);
            prodsGs_[iAux]->push_back(prodG);
        }
        clog << "Initialized prodsGs[" << iAux << "].\n";
        printVecVec(prodsGs_, "prodsGs");
    }

    /*! Initializes obstacle set for all layers of the base system in baseOs_.
     *  \param[in] addO     Pointer to function that adds obstacle states to a SymbolicSet.
     */
    template<class O_t>
    void initializeBaseOs(O_t addO) {
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicSet* baseO = new SymbolicSet(*baseXs_[iAbs]);
            addO(baseO);
            baseOs_.push_back(baseO);
        }
        clog << "Initialized baseOs.\n";
    }

    /*! Takes base system and auxiliary systems' transition systems and yields product transition systems for each layer. */
    void composeAbstractions() {
        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                ((*prodsTs_[iAux])[iAbs])->symbolicSet_ = baseTs_[iAbs]->symbolicSet_ & ((*auxsTs_[iAux])[iAbs])->symbolicSet_;
                ((*prodsTTs_[iAux])[iAbs])->symbolicSet_ = ((*prodsTs_[iAux])[iAbs])->symbolicSet_.ExistAbstract(*((*prodsNotXUVars_[iAux])[iAbs]));
            }
        }
        clog << "Composed transition relations, stored in prods' Ts.\n";
        clog << "Composed reduced transition relations, stored in prods' TTs.\n";

        printVecVec(prodsTs_, "prodsTs");
        printVecVec(prodsTTs_, "prodsTTs");
    }

    /*! Initializes data members. */
    void initialize(System* base, vector<System*> auxs) {
        ddmgr_ = new Cudd;
        base_ = base;
        auxs_ = auxs;

        initializeVecVecs();
        initializeEtaTau();
        initializeSolvers();
        initializeBDDs();
    }

    /*! Computes transition systems for the base system in each layer by iteratively sampling over all state-input pairs.
     *  \param[in]  sysNext     Function pointer that, given a state and an input, modifies that state according to system dynamics.
     *  \param[in]  radNext     Function pointer that, given the initial growth bound (eta/2), modifies it according to the growth bound equation.
     *  \param[in]  x           Only used for passing type of a state.
     *  \param[in]  u           Only used for passing type of an input.
     */
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
        clog << "Computed base's transition relations.\n";
    }

    /*! Computes transition systems for auxiliary systems in each layer by iteratively sampling over all states. A dummy input space containing one element is used for the computation, then abstracted out.
     *  \param[in]  sysNext     Function pointer that, given a state (and a dummy input), modifies that state according to system dynamics.
     *  \param[in]  radNext     Function pointer that, given the initial growth bound (eta/2), modifies it according to the growth bound equation.
     *  \param[in]  x           Only used for passing type of a state.
     *  \param[in]  u           Only used for passing type of a dummy input.
     */
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
        clog << "Computed aux" << iAux << "'s transition relations.\n";
    }

    /*! Initializes Runge-Kutta numerical approximators for use in computing the abstract transition systems. */
    void initializeSolvers() {
        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            vector<OdeSolver*>* auxSolvers = new vector<OdeSolver*>;
            auxsSolvers_.push_back(auxSolvers);
        }

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            OdeSolver* baseSolver = new OdeSolver(*base_->dimX_, *base_->nSubInt_, taus_[iAbs][0]);
            baseSolvers_.push_back(baseSolver);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                System* aux = auxs_[iAux];
                OdeSolver* auxSolver = new OdeSolver(*aux->dimX_, *base_->nSubInt_, taus_[iAbs][0]);
                auxsSolvers_[iAux]->push_back(auxSolver);
            }
        }

        clog << "Initialized base's, auxs' ODE solvers.\n";
    }

    /*! Initializes grid parameters and time sampling parameters for each layer. */
    void initializeEtaTau() {

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
            for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
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

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            vector<double*>* prodEtaX = new vector<double*>;
            for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
                double* prodEtai = new double[*base_->dimX_ + *auxs_[iAux]->dimX_];
                for (int iDim = 0; iDim < *base_->dimX_; iDim++) {
                    prodEtai[iDim] = baseEtaXs_[iAbs][iDim];
                }
                for (int iDim = 0; iDim < *auxs_[iAux]->dimX_; iDim++) {
                    prodEtai[iDim + *base_->dimX_] = (*(auxsEtaXs_[iAux]))[iAbs][iDim];
                }
                prodEtaX->push_back(prodEtai);
            }
            prodsEtaXs_.push_back(prodEtaX);
        }

        double* allTauCur = new double;
        *allTauCur = *base_->tau_;

        for (int i = 0; i < *base_->numAbs_; i++) {
            double* allTaui = new double;
            *allTaui = *allTauCur;
            taus_.push_back(allTaui);
            *allTauCur /= *base_->tauRatio_;
        }
        delete allTauCur;

        clog << "Initialized base's, auxs', prods' etaXs, taus.\n";

        printVecArray(baseEtaXs_, "baseEtaXs_", *base_->dimX_);
        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            string prefix = "auxsEtaXs_[";
            prefix += std::to_string(iAux);
            prefix += "]";
            printVecArray(*auxsEtaXs_[iAux], prefix, *auxs_[iAux]->dimX_);
        }

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            string prefix = "prodsEtaXs_[";
            prefix += std::to_string(iAux);
            prefix += "]";
            printVecArray(*prodsEtaXs_[iAux], prefix, *base_->dimX_ + *auxs_[iAux]->dimX_);
        }

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            double* prodEtaRatio = new double[*base_->dimX_ + *auxs_[iAux]->dimX_];
            for (int iDim = 0; iDim < *base_->dimX_; iDim++) {
                prodEtaRatio[iDim] = base_->etaRatio_[iDim];
            }
            for (int iDim = 0; iDim < *auxs_[iAux]->dimX_; iDim++) {
                prodEtaRatio[iDim + *base_->dimX_] = auxs_[iAux]->etaRatio_[iDim];
            }
            prodsEtaRatio_.push_back(prodEtaRatio);

            clog << "prodsEtaRatio[" << iAux << "]: ";
            printArray(prodEtaRatio, *base_->dimX_ + *auxs_[iAux]->dimX_);
        }
        clog << "Initialized prods' etaRatio.\n";
    }

    /*! Initializes vectors of vectors of SymbolicSets, BDDs, and others. Each vector of a SymbolicSet or BDD is always for a particular product system.
     *  For all such vectors excluding vectors of prodsFinalZs_ and prodsFinalCs_, each element is for a particular abstraction layer, from coarsest to finest.
     */
    void initializeVecVecs() {
        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            vector<SymbolicSet*>* auxXs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* auxX2s = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* auxXXs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodXs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodX2s = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodXXs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodZs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodYs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodPreYs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodPrevYs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodPrevPreYs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodCs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodTs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodTTs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodGs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* auxTs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodValidZs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodValidCs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodFinalZs = new vector<SymbolicSet*>;
            vector<SymbolicSet*>* prodFinalCs = new vector<SymbolicSet*>;
            auxsXs_.push_back(auxXs);
            auxsX2s_.push_back(auxX2s);
            auxsXXs_.push_back(auxXXs);
            prodsXs_.push_back(prodXs);
            prodsX2s_.push_back(prodX2s);
            prodsXXs_.push_back(prodXXs);
            prodsZs_.push_back(prodZs);
            prodsYs_.push_back(prodYs);
            prodsPreYs_.push_back(prodPreYs);
            prodsPrevYs_.push_back(prodPrevYs);
            prodsPrevPreYs_.push_back(prodPrevPreYs);
            prodsCs_.push_back(prodCs);
            prodsTs_.push_back(prodTs);
            prodsTTs_.push_back(prodTTs);
            prodsGs_.push_back(prodGs);
            auxsTs_.push_back(auxTs);
            prodsValidZs_.push_back(prodValidZs);
            prodsValidCs_.push_back(prodValidCs);
            prodsFinalZs_.push_back(prodFinalZs);
            prodsFinalCs_.push_back(prodFinalCs);

            vector<BDD*>* auxCubesX = new vector<BDD*>;
            vector<BDD*>* auxCubesX2 = new vector<BDD*>;
            vector<int*>* prodPermutesXtoX2 = new vector<int*>;
            vector<int*>* prodPermutesX2toX = new vector<int*>;
            vector<BDD*>* prodNotXUVars = new vector<BDD*>;
            vector<BDD*>* prodNotXVars = new vector<BDD*>;
            vector<int*>* prodPermutesCoarser = new vector<int*>;
            vector<int*>* prodPermutesFiner = new vector<int*>;
            vector<BDD*>* prodCubeCoarser = new vector<BDD*>;
            auxsCubesX_.push_back(auxCubesX);
            auxsCubesX2_.push_back(auxCubesX2);
            prodsPermutesXtoX2_.push_back(prodPermutesXtoX2);
            prodsPermutesX2toX_.push_back(prodPermutesX2toX);
            prodsNotXUVars_.push_back(prodNotXUVars);
            prodsNotXVars_.push_back(prodNotXVars);
            prodsPermutesCoarser_.push_back(prodPermutesCoarser);
            prodsPermutesFiner_.push_back(prodPermutesFiner);
            prodsCubesCoarser_.push_back(prodCubeCoarser);
        }
    }

    /*! Initializes SymbolicSets, BDDs, and others related to symbolic representation of abstractions. */
    void initializeBDDs() {

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            SymbolicSet* baseX = new SymbolicSet(*ddmgr_, *base_->dimX_, base_->lbX_, base_->ubX_, baseEtaXs_[iAbs], taus_[iAbs][0]);
            baseX->addGridPoints();
            baseXs_.push_back(baseX);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                System* aux = auxs_[iAux];
                SymbolicSet* auxX = new SymbolicSet(*ddmgr_, *aux->dimX_, aux->lbX_, aux->ubX_, (*auxsEtaXs_[iAux])[iAbs], taus_[iAbs][0]);
                auxX->addGridPoints();
                auxsXs_[iAux]->push_back(auxX);

                SymbolicSet* prodX = new SymbolicSet(*baseXs_[iAbs], *auxX);
                prodX->addGridPoints();
                prodsXs_[iAux]->push_back(prodX);

                SymbolicSet* prodZ = new SymbolicSet(*prodX);
                prodsZs_[iAux]->push_back(prodZ);

                SymbolicSet* prodY = new SymbolicSet(*prodX);
                prodsYs_[iAux]->push_back(prodY);

                SymbolicSet* prodPreY = new SymbolicSet(*prodX);
                prodsPreYs_[iAux]->push_back(prodPreY);

                SymbolicSet* prodPrevY = new SymbolicSet(*prodX);
                prodsPrevYs_[iAux]->push_back(prodPrevY);

                SymbolicSet* prodPrevPreY = new SymbolicSet(*prodX);
                prodsPrevPreYs_[iAux]->push_back(prodPrevPreY);

                SymbolicSet* prodValidZ = new SymbolicSet(*prodX);
                prodsValidZs_[iAux]->push_back(prodValidZ);
            }
        }
        clog << "Initialized base's Xs to full.\n";
        clog << "Initialized auxs' Xs to full.\n";
        clog << "Initialized prods' Xs to full, Ys to empty, preYs to empty, prevYs to empty, prevPreYs to empty, Zs to empty, validZs to empty.\n";

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
        clog << "Initialized base's X2s to empty.\n";
        clog << "Initialized auxs' X2s to empty.\n";
        clog << "Initialized prods' X2s to empty.\n";


        for (int iAbs = 0; iAbs < *base_->numAbs_ - 1; iAbs++) {
            SymbolicSet* baseXX = new SymbolicSet(*baseXs_[iAbs], *baseXs_[iAbs+1]);
            mapAbstractions(baseXs_[iAbs], baseXs_[iAbs+1], baseXX);
            baseXXs_.push_back(baseXX);

            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* auxXX = new SymbolicSet(*((*auxsXs_[iAux])[iAbs]), *((*auxsXs_[iAux])[iAbs+1]));
                mapAbstractions((*auxsXs_[iAux])[iAbs], (*auxsXs_[iAux])[iAbs+1], auxXX);
                auxsXXs_[iAux]->push_back(auxXX);

                SymbolicSet* prodXX = new SymbolicSet(*((*prodsXs_[iAux])[iAbs]), *((*prodsXs_[iAux])[iAbs+1]));
                prodXX->symbolicSet_ = baseXX->symbolicSet_ & auxXX->symbolicSet_;
                prodsXXs_[iAux]->push_back(prodXX);
            }
        }
        clog << "Initialized base's, auxs', via mapAbstractions, and prods' XXs via composition.\n";

        baseU_ = new SymbolicSet(*ddmgr_, *base_->dimU_, base_->lbU_, base_->ubU_, base_->etaU_, 0);
        baseU_->addGridPoints();

        System* aux = auxs_[0];
        auxU_ = new SymbolicSet(*ddmgr_, *aux->dimU_, aux->lbU_, aux->ubU_, aux->etaU_, 0);
        auxU_->addGridPoints();

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* prodC = new SymbolicSet(*((*prodsXs_[iAux])[iAbs]), *baseU_);
                prodsCs_[iAux]->push_back(prodC);

                SymbolicSet* prodValidC = new SymbolicSet(*prodC);
                prodsValidCs_[iAux]->push_back(prodValidC);
            }
        }
        clog << "Initialized prods' Cs to empty, validCs to empty.\n";

        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                SymbolicSet* prodT = new SymbolicSet(*((*prodsCs_[iAux])[iAbs]), *((*prodsX2s_[iAux])[iAbs]));
                prodsTs_[iAux]->push_back(prodT);
                SymbolicSet* prodTT = new SymbolicSet(*((*prodsCs_[iAux])[iAbs]));
                prodsTTs_[iAux]->push_back(prodTT);
            }
        }
        clog << "Initialized prods' Ts, TTs to empty.\n";

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

                *auxCubeX = ((*auxsXs_[iAux])[iAbs])->getCube();
                *auxCubeX2 = ((*auxsX2s_[iAux])[iAbs])->getCube();

                auxsCubesX_[iAux]->push_back(auxCubeX);
                auxsCubesX2_[iAux]->push_back(auxCubeX2);
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
        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            string prefix = "prodsPermutesXtoX2_[";
            prefix += std::to_string(iAux);
            prefix += "]";
            printVecArray(*prodsPermutesXtoX2_[iAux], prefix, numBDDVars);
        }


        for (int iAbs = 0; iAbs < *base_->numAbs_; iAbs++) {
            for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
                BDD* prodNotXUVar = new BDD;
                BDD* prodNotXVar = new BDD;
                *prodNotXUVar = ddmgr_->bddOne();

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
            *baseNotXUVar = ddmgr_->bddOne();

            *baseNotXUVar &= *((*prodsNotXUVars_[0])[iAbs]);
            *baseNotXUVar &= *((*auxsCubesX_[0])[iAbs]);
            *baseNotXVar = *baseNotXUVar & baseU_->getCube();
            baseNotXUVars_.push_back(baseNotXUVar);
            baseNotXVars_.push_back(baseNotXVar);
        }
        clog << "Initialized BDDs used for projection (ExistAbstract).\n";

        int* ones = new int[auxs_.size()]; // number of ones in each prodXs' etaRatio

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            for (size_t iDim = 0; iDim < (*prodsXs_[iAux])[0]->dim_; iDim++) {
                if (prodsEtaRatio_[iAux][iDim] == 1) {
                    ones[iAux] = ones[iAux] + 1;
                }
            }
        }
        clog << "here\n";

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            for (int iAbs = 1; iAbs < *base_->numAbs_; iAbs++) {
                BDD* varsCoarser = new BDD[(*prodsXs_[iAux])[iAbs]->dim_ - ones[iAux]];
                int ind = 0;
                for (size_t iDim = 0; iDim < (*prodsXs_[iAux])[iAbs]->dim_; iDim++) {
                    if (prodsEtaRatio_[iAux][iDim] == 2) {
                        varsCoarser[ind] = ddmgr_->bddVar((*prodsXs_[iAux])[iAbs]->indBddVars_[iDim][0]);
                        ind++;
                    }
                }
                BDD* cubeCoarser = new BDD;
                *cubeCoarser = ddmgr_->bddComputeCube(varsCoarser, NULL, (*prodsXs_[iAux])[iAbs]->dim_ - ones[iAux]);
                prodsCubesCoarser_[iAux]->push_back(cubeCoarser);
                delete[] varsCoarser;
            }
        }
        clog << "Initialized prodsCubesCoarser.\n";

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            for (int iAbs = 1; iAbs < *base_->numAbs_; iAbs++) {
                int* permuteCoarser = new int[numBDDVars];
                for (int ind = 0; ind < numBDDVars; ind++) {
                    permuteCoarser[ind] = 0;
                }
                for (size_t iDim = 0; iDim < (*prodsXs_[iAux])[iAbs]->dim_; iDim++) {
                    if (prodsEtaRatio_[iAux][iDim] == 2) {
                        for (size_t projInd = 1; projInd < (*prodsXs_[iAux])[iAbs]->nofBddVars_[iDim]; projInd++) {
                            permuteCoarser[(*prodsXs_[iAux])[iAbs]->indBddVars_[iDim][projInd]] = (*prodsXs_[iAux])[iAbs-1]->indBddVars_[iDim][projInd-1];
                        }
                    }
                    else if (prodsEtaRatio_[iAux][iDim] == 1) {
                        for (size_t projInd = 0; projInd < (*prodsXs_[iAux])[iAbs]->nofBddVars_[iDim]; projInd++) {
                            permuteCoarser[(*prodsXs_[iAux])[iAbs]->indBddVars_[iDim][projInd]] = (*prodsXs_[iAux])[iAbs-1]->indBddVars_[iDim][projInd];
                        }
                    }
                }
                prodsPermutesCoarser_[iAux]->push_back(permuteCoarser);

                clog << "permuteCoarser[" << iAux << "][" << iAbs << "]: ";
                printArray(permuteCoarser, numBDDVars);
            }
        }

        for (size_t iAux = 0; iAux < auxs_.size(); iAux++) {
            for (int iAbs = 0; iAbs < *base_->numAbs_ - 1; iAbs++) {
                int* permuteFiner = new int[numBDDVars];
                for (int ind = 0; ind < numBDDVars; ind++) {
                    permuteFiner[ind] = 0;
                }

                for (size_t iDim = 0; iDim < (*prodsXs_[iAux])[iAbs]->dim_; iDim++) {
                    if (prodsEtaRatio_[iAux][iDim] == 2) {
                        for (size_t projInd = 0; projInd < (*prodsXs_[iAux])[iAbs]->nofBddVars_[iDim]; projInd++) {
                            permuteFiner[(*prodsXs_[iAux])[iAbs]->indBddVars_[iDim][projInd]] = (*prodsXs_[iAux])[iAbs+1]->indBddVars_[iDim][projInd+1];
                        }
                    }
                    else if (prodsEtaRatio_[iAux][iDim] == 1) {
                        for (size_t projInd = 0; projInd < (*prodsXs_[iAux])[iAbs]->nofBddVars_[iDim]; projInd++) {
                            permuteFiner[(*prodsXs_[iAux])[iAbs]->indBddVars_[iDim][projInd]] = (*prodsXs_[iAux])[iAbs+1]->indBddVars_[iDim][projInd];
                        }
                    }
                }
                prodsPermutesFiner_[iAux]->push_back(permuteFiner);

                clog << "permuteFiner[" << iAux << "][" << iAbs << "]: ";
                printArray(permuteFiner, numBDDVars);
            }
        }

        delete[] ones;



        printVec(baseXs_, "baseXs");
        //        printVec(baseXXs_, "baseXXs");

        //        printVecVec(auxsXs_, "auxsXs");
        //        printVecVec(auxsXXs_, "auxsXXs");

        printVecVec(prodsXs_, "prodsXs");
        printVecVec(prodsXXs_, "prodsXXs");

        //        printVecVec(prodsTs_, "prodsTs");


//        checkNotVars();
    }

    /*! Generates a mapping between two consecutive state space abstractions.
     *  \param[in]      Xc          Coarser state space abstraction.
     *  \param[in]      Xf          Finer state space abstraction.
     *  \param[in,out]  XX    Mapping for a finer cell to the coarser cell that is its superset.
     */
    void mapAbstractions(SymbolicSet* Xc, SymbolicSet* Xf, SymbolicSet* XX) {


        int* XfMinterm;
        double* xPoint = new double[Xc->dim_];
        vector<double> XcPoint (Xc->dim_, 0);
        double* XXPoint = new double[XX->dim_];

        int totalIter = Xf->symbolicSet_.CountMinterm(Xf->nvars_);
        int iter = 0;

        cout << "XXs totalIter: " << totalIter << '\n';

        for (Xf->begin(); !Xf->done(); Xf->next()) {
            iter++;
            if (iter % 50000 == 0) {
                cout << iter << "\n";
            }

            XfMinterm = (int*)Xf->currentMinterm();
            Xf->mintermToElement(XfMinterm, xPoint);
            for (size_t i = 0; i < Xc->dim_; i++) {
                XcPoint[i] = xPoint[i];
            }
            if (!(Xc->isElement(XcPoint))) {
                continue;
            }
            for (size_t i = 0; i < XX->dim_; i++) {
                XXPoint[i] = xPoint[i % Xc->dim_];
            }
            XX->addPoint(XXPoint);
        }
        cout << '\n';

        delete[] xPoint;
        delete[] XXPoint;
    }

    /*! Debugging function. */
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

#endif /* Composition_HH_ */
