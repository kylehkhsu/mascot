/*! \file GBuchi.hh
 *  Contains the GBuchi class. */

#ifndef GBUCHI_HH_
#define GBUCHI_HH_

#include "Composition.hh"

using std::vector;
using std::cout;

namespace scots {

/*! \class GBuchi
 *  \brief A class (derived from base Composition) that does adaptive multiscale abstraction-based synthesis for a generalized Buchi specification.
 */
class GBuchi: public Composition {
public:
    int startAbs_;
    int minToGoCoarser_;
    int minToBeValid_;
    int earlyBreak_;
    int verbose_;



    /*! Constructor for a GBuchi object. */
    GBuchi(char* logFile)
        : Composition(logFile)
    {
    }

//    template<class I_type>
//    void auxAddInitial(I_type addI) {

//    }

    ~GBuchi() {

    }

    void gBuchi(int startAbs, int minToGoCoarser, int minToBeValid, int verbose = 1) {
        startAbs_ = startAbs;
        minToGoCoarser_ = minToGoCoarser;
        minToBeValid_ = minToBeValid;
        earlyBreak_ = 0;
        verbose_ = verbose;

        int fAbs = *this->base_->numAbs_ - 1;

        if ( (startAbs_ < 0) || (startAbs_ >= *this->base_->numAbs_) ) {
            error("Error: startAbs out of range.\n");
        }

        TicToc tt;
        tt.tic();

        // initialization
        size_t prevAux = this->auxs_.size() - 1;
        size_t curAux = 0;
        size_t numConv = 0;
        int nuIter = 1;
        int nuStop = 0;

        for (size_t iAux = 0; iAux < this->auxs_.size(); iAux++) {
            ((*this->prodsYs_[iAux])[fAbs])->addGridPoints(); // Y = Q
            ((*this->prodsPrevYs_[iAux])[fAbs])->symbolicSet_ = this->ddmgr_->bddZero();
        }

        while (1) {
            nuMu(&prevAux, &curAux, &numConv, &nuIter, &nuStop);
            if (nuStop) {
                break;
            }
        }

        clog << "----------------------------------------generalizedBuchi: ";
        tt.toc();

        checkMakeDir("C");
        checkMakeDir("Z");
        checkMakeDir("G");
        string prefix = "";
        for (size_t iAux = 0; iAux < this->auxs_.size(); iAux++) {
            prefix = "C/C";
            prefix += std::to_string(iAux+1);
            saveVec(*prodsFinalCs_[iAux], prefix);
            prefix = "Z/Z";
            prefix += std::to_string(iAux+1);
            saveVec(*prodsFinalZs_[iAux], prefix);
            prefix = "G/G";
            prefix += std::to_string(iAux+1);
            saveVec(*prodsGs_[iAux], prefix);
        }

        checkMakeDir("plotting");
        this->baseXs_[0]->writeToFile("plotting/X.bdd");
        this->baseOs_[fAbs]->writeToFile("plotting/O.bdd");

    }

    void alwaysEventuallyOne(int startAbs, int minToGoCoarser, int minToBeValid, int verbose = 1) {
        startAbs_ = startAbs;
        minToGoCoarser_ = minToGoCoarser;
        minToBeValid_ = minToBeValid;
        earlyBreak_ = 0;
        verbose_ = verbose;

        if ( (startAbs_ < 0) || (startAbs_ >= *this->base_->numAbs_) ) {
            error("Error: startAbs out of range.\n");
        }
        if (this->auxs_.size() != 1) {
            error("Error: this debugging function only works with one auxiliary system.");
        }

        TicToc tt;
        tt.tic();

        size_t prevAux = 0;
        size_t curAux = 0;
        size_t numConv = 0;
        int nuIter = 1;
        int nuStop = 0;

        // initialization for nuMu
        ((*this->prodsYs_[curAux])[*this->base_->numAbs_ - 1])->addGridPoints(); // Y = Q
        ((*this->prodsPrevYs_[curAux])[*this->base_->numAbs_ - 1])->symbolicSet_ = this->ddmgr_->bddZero();


        while (1) {
            nuMu(&prevAux, &curAux, &numConv, &nuIter, &nuStop);
            if (nuStop) {
                break;
            }
        }

        int fAbs = *this->base_->numAbs_ - 1;

        SymbolicSet* lastZ = (this->prodsFinalZs_[curAux])->back();

        // part of goal that we can reach infinitely often
        ((*this->prodsGs_[curAux])[fAbs])->symbolicSet_ &= lastZ->symbolicSet_;

        checkMakeDir("C");
        string prefix = "C/C";
        prefix += std::to_string(curAux+1);
        saveVec(*prodsFinalCs_[curAux], prefix);
        checkMakeDir("Z");
        prefix = "Z/Z";
        prefix += std::to_string(curAux+1);
        saveVec(*prodsFinalZs_[curAux], prefix);
        checkMakeDir("plotting");
        this->baseXs_[0]->writeToFile("plotting/X.bdd");
        this->baseOs_[fAbs]->writeToFile("plotting/O.bdd");
        ((*this->prodsIs_[curAux])[fAbs])->writeToFile("plotting/I.bdd");
        checkMakeDir("G");
        ((*this->prodsGs_[curAux])[fAbs])->writeToFile("G/G.bdd");


        clog << "----------------------------------------alwaysEventually: ";
        tt.toc();
    }




    void nuMu(size_t* prevAux, size_t* curAux, size_t* numConv, int* nuIter, int* nuStop) {
        clog << "------------------------------nuMu iteration: " << *nuIter << "------------------------------\n";
        cout << "------------------------------nuMu iteration: " << *nuIter << "------------------------------\n";

        // find pre(Y) of product system prevAux at finest abstraction
        int fAbs = *this->base_->numAbs_ - 1;
        BDD preYC = Cpre(((*this->prodsYs_[*prevAux])[fAbs])->symbolicSet_,
                         ((*this->prodsTs_[*prevAux])[fAbs])->symbolicSet_,
                         ((*this->prodsTTs_[*prevAux])[fAbs])->symbolicSet_,
                         ((*this->prodsPermutesXtoX2_[*prevAux])[fAbs]),
                         *((*this->prodsNotXUVars_[*prevAux])[fAbs]));
        BDD basePreY = preYC.ExistAbstract(*this->baseNotXVars_[fAbs]);

        // obtain pre(Y) to use in this mu.
        ((*this->prodsPreYs_[*curAux])[fAbs])->symbolicSet_ = basePreY & ((*this->auxsXs_[*curAux])[fAbs])->symbolicSet_;

        cout << "prodPreY finest:\n";
        ((*this->prodsPreYs_[*curAux])[fAbs])->printInfo(1);

        // construct pre(Y) for all other abstractions
        for (int iAbs = fAbs; iAbs > 0; iAbs--) {
            // reset to 0 since innerCoarser only adds states
            ((*this->prodsPreYs_[*curAux])[iAbs-1])->symbolicSet_ = this->ddmgr_->bddZero();
            innerCoarser((*this->prodsPreYs_[*curAux])[iAbs-1],
                    (*this->prodsPreYs_[*curAux])[iAbs],
                    (*this->prodsXXs_[*curAux])[iAbs-1],
                    (*this->prodsNotXVars_[*curAux])[iAbs-1],
                    (*this->prodsXs_[*curAux])[iAbs-1],
                    this->prodsNumFiner_[*curAux]);
        }

        // initialization for mu
        int curAbs = startAbs_;
        int muIter = 1;
        int justCoarsed = 0;
        int iterCurAbs = 1;
        int reached = 0;
        int muStop = 0;

        for (int iAbs = 0; iAbs < *this->base_->numAbs_; iAbs++) {
            ((*this->prodsZs_[*curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
            ((*this->prodsValidZs_[*curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
            ((*this->prodsCs_[*curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
            ((*this->prodsValidCs_[*curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
        }

        clog << "Resetting prodCs, prodZs, prodValidCs, prodValidZs to empty.\n";
        clog << "Emptying prodFinalCs, prodFinalZs.\n";
        size_t j = prodsFinalCs_[*curAux]->size();
        for (size_t i = 0; i < j; i++) {
            delete prodsFinalCs_[*curAux]->back();
            prodsFinalCs_[*curAux]->pop_back();
            delete prodsFinalZs_[*curAux]->back();
            prodsFinalZs_[*curAux]->pop_back();
        }

        while (1) { // eventually
            this->mu(*curAux, &curAbs, &muIter, &justCoarsed, &iterCurAbs, &reached, &muStop);
            if (muStop) {
                break;
            }
        }

        // projecting all winning states onto finest abstraction
        ((*this->prodsYs_[*curAux])[0])->symbolicSet_ = ((*this->prodsZs_[*curAux])[0])->symbolicSet_;
        for (int iAbs = 0; iAbs < *this->base_->numAbs_ - 1; iAbs++) {
            innerFiner((*this->prodsYs_[*curAux])[iAbs],
                    (*this->prodsYs_[*curAux])[iAbs+1],
                    (*this->prodsXXs_[*curAux])[iAbs],
                    (*this->prodsNotXVars_[*curAux])[iAbs+1],
                    (*this->prodsXs_[*curAux])[iAbs+1]);
            ((*this->prodsYs_[*curAux])[iAbs+1])->symbolicSet_ |= ((*this->prodsZs_[*curAux])[iAbs+1])->symbolicSet_;
        }

        if ( ((*this->prodsYs_[*curAux])[fAbs])->symbolicSet_ == ((*this->prodsPrevYs_[*curAux])[fAbs])->symbolicSet_ ) {
            *numConv += 1;
        }
        else {
            *numConv = 0;
            ((*this->prodsPrevYs_[*curAux])[fAbs])->symbolicSet_ = ((*this->prodsYs_[*curAux])[fAbs])->symbolicSet_;
        }
        if (*numConv == this->auxs_.size()) {
            *nuStop = 1;
            return;
        }

        *nuIter += 1;

        *prevAux = *curAux;
        if (*curAux == this->auxs_.size() - 1) {
            *curAux = 0;
        }
        else {
            *curAux += 1;
        }
    }



    void mu(int curAux, int* curAbs, int* iter, int* justCoarsed, int* iterCurAbs, int* reached, int* stop) {
        clog << "current abstraction: " << *curAbs << '\n';
        clog << "mu iteration: " << *iter << '\n';
        clog << "justCoarsed: " << *justCoarsed << '\n';
        clog << "iterCurAbs: " << *iterCurAbs << '\n';
        clog << "reached: " << *reached << '\n';
        clog << "controllers: " << this->prodsFinalCs_[curAux]->size() << '\n';
        cout << "current abstraction: " << *curAbs << '\n';
        cout << "mu iteration: " << *iter << '\n';
        cout << "justCoarsed: " << *justCoarsed << '\n';
        cout << "iterCurAbs: " << *iterCurAbs << '\n';
        cout << "reached: " << *reached << '\n';
        cout << "controllers: " << this->prodsFinalCs_[curAux]->size() << '\n';

        // pre(Z)
        BDD C = Cpre(((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_,
                ((*this->prodsTs_[curAux])[*curAbs])->symbolicSet_,
                ((*this->prodsTTs_[curAux])[*curAbs])->symbolicSet_,
                ((*this->prodsPermutesXtoX2_[curAux])[*curAbs]),
                *((*this->prodsNotXUVars_[curAux])[*curAbs]));
        // pre(Z) \cup (G \cap pre(Y))
        C |= ( (((*this->prodsGs_[curAux])[*curAbs])->symbolicSet_) & (((*this->prodsPreYs_[curAux])[*curAbs])->symbolicSet_) );
        // new {(x,u)}
        BDD N = C & (!((*this->prodsCs_[curAux])[*curAbs])->symbolicSet_.ExistAbstract(*((*this->prodsNotXVars_[curAux])[*curAbs])));
        // add new {(x,u)} to C
        ((*this->prodsCs_[curAux])[*curAbs])->symbolicSet_ |= N;
        // project C onto Z
        ((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_ = C.ExistAbstract(*((*this->prodsNotXVars_[curAux])[*curAbs]));

        if ( (((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_ & ((*this->prodsIs_[curAux])[*curAbs])->symbolicSet_) != this->ddmgr_->bddZero() ) {
            *reached = 1;
            if (earlyBreak_ == 1) {
                saveCZ(curAux, *curAbs);
                *stop = 1;
                clog << "\nmu number of controllers: " << this->prodsFinalCs_[curAux]->size() << '\n';
                return;
            }
        }

        if (*iter != 1) {
            if (N == this->ddmgr_->bddZero()) {
                if (*curAbs == *this->base_->numAbs_ - 1) {
                    saveCZ(curAux, *curAbs);
                    if (*reached == 1) {
                        clog << "Reached.\n";
                    }
                    else {
                        clog << "Did not reach.\n";
                    }
                    *stop = 1;
                    clog << "\nmu number of controllers: " << this->prodsFinalCs_[curAux]->size() << '\n';
                }
                else {
                    if (verbose_) {
                        clog << "No new winning states; going finer.\n";
                    }
                    if (*justCoarsed == 1) {
                        if (verbose_) {
                            clog << "Current controller has not been declared valid.\n";
                            clog << "Removing last elements of prodsFinalCs_[" << curAux << "] and prodsFinalZs_[" << curAux << "].\n";
                            clog << "Resetting prodsZs_[" << curAux << "][" << *curAbs << "] and prodsCs_[" << curAux << "][" << *curAbs << "] from valids.\n";
                        }

                        SymbolicSet* prodC = prodsFinalCs_[curAux]->back();
                        SymbolicSet* prodZ = prodsFinalZs_[curAux]->back();
                        prodsFinalCs_[curAux]->pop_back();
                        prodsFinalZs_[curAux]->pop_back();
                        delete(prodC);
                        delete(prodZ);

                        ((*this->prodsCs_[curAux])[*curAbs])->symbolicSet_ = ((*this->prodsValidCs_[curAux])[*curAbs])->symbolicSet_;
                        ((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_ = ((*this->prodsValidZs_[curAux])[*curAbs])->symbolicSet_;
                        *justCoarsed = 0;
                    }
                    else {
                        if (verbose_) {
                            clog << "Current controller has been declared valid.\n";
                            clog << "Saving prodsZs_[" << curAux << "][" << *curAbs << "] and prodsCs_[" << curAux << "][" << *curAbs << "] into prodValids, prodFinals.\n";
                            clog << "Finding inner approximation of prodsZs_[" << curAux << "][" << *curAbs << "] and copying into prodsZs_[" << curAux << "][" << *curAbs+1 << "].\n";
                        }
                        ((*this->prodsValidCs_[curAux])[*curAbs])->symbolicSet_ = ((*this->prodsCs_[curAux])[*curAbs])->symbolicSet_;
                        ((*this->prodsValidZs_[curAux])[*curAbs])->symbolicSet_ = ((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_;

                        saveCZ(curAux, *curAbs);

                        innerFiner((*this->prodsZs_[curAux])[*curAbs],
                                (*this->prodsZs_[curAux])[*curAbs+1],
                                (*this->prodsXXs_[curAux])[*curAbs],
                                (*this->prodsNotXVars_[curAux])[*curAbs+1],
                                (*this->prodsXs_[curAux])[*curAbs+1]);
                        ((*this->prodsValidZs_[curAux])[*curAbs+1])->symbolicSet_ = ((*this->prodsZs_[curAux])[*curAbs+1])->symbolicSet_;
                    }
                    *curAbs += 1;
                    *iterCurAbs = 1;
                }
            }
            else {
                if ((*justCoarsed == 1) && (*iterCurAbs >= minToBeValid_)) {
                    if (verbose_) {
                        clog << "Current controller now valid.\n";
                    }
                    *justCoarsed = 0;
                }
                if (*curAbs == 0) {
                    if (verbose_) {
                        clog << "More new states, already at coarsest gridding, continuing.\n";
                    }
                    *iterCurAbs += 1;
                }
                else {
                    if (*iterCurAbs >= minToGoCoarser_) {
                        if (verbose_) {
                            clog << "More new states, minToGoCoarser achieved; try going coarser.\n";
                        }
                        int more = innerCoarser((*this->prodsZs_[curAux])[*curAbs-1],
                                (*this->prodsZs_[curAux])[*curAbs],
                                (*this->prodsXXs_[curAux])[*curAbs-1],
                                (*this->prodsNotXVars_[curAux])[*curAbs-1],
                                (*this->prodsXs_[curAux])[*curAbs-1],
                                this->prodsNumFiner_[curAux]);
                        if (more == 0) {
                            if (verbose_) {
                                clog << "Projecting to coarser gives no more states in coarser; not changing abstraction.\n";
                            }
                            *iterCurAbs += 1;
                        }
                        else {
                            if (verbose_) {
                                clog << "Projecting to coarser gives more states in coarser.\n";
                                clog << "Saving prodsZs_[" << curAux << "][" << *curAbs << "] and prodsCs_[" << curAux << "][" << *curAbs << "] into prodValids, prodFinals.\n";
                                clog << "Saving prodsZs_[" << curAux << "][" << *curAbs-1 << "] and prodsCs_[" << curAux << "][" << *curAbs-1 << "] into prodValids.\n";
                            }
                            saveCZ(curAux, *curAbs);
                            ((*this->prodsValidCs_[curAux])[*curAbs])->symbolicSet_ = ((*this->prodsCs_[curAux])[*curAbs])->symbolicSet_;
                            ((*this->prodsValidZs_[curAux])[*curAbs])->symbolicSet_ = ((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_;

                            ((*this->prodsValidCs_[curAux])[*curAbs-1])->symbolicSet_ = ((*this->prodsCs_[curAux])[*curAbs-1])->symbolicSet_;
                            ((*this->prodsValidZs_[curAux])[*curAbs-1])->symbolicSet_ = ((*this->prodsZs_[curAux])[*curAbs-1])->symbolicSet_;

                            *justCoarsed = 1;
                            *curAbs -= 1;
                            *iterCurAbs = 1;
                        }
                    }
                    else {
                        *iterCurAbs += 1;
                    }
                }
            }
        }
        else {
            if (verbose_) {
                clog << "First iteration, nothing happens.\n";
            }
        }
        clog << "\n";
        *iter += 1;
    }

    BDD Cpre (BDD Z, BDD T, BDD TT, int* permuteXtoX2, BDD notXUVar) {
        // swap to post-state variables
        BDD Z2 = Z.Permute(permuteXtoX2);
        // non-winning post-states
        BDD nZ2 = !Z2;
        // {(x,u)} with non-winning post-state
        BDD nC = T.AndAbstract(nZ2, notXUVar);
        // {(x,u)} with no non-winning post-states
        BDD C = !nC;
        // get rid of non-domain information
        C = TT.AndAbstract(C, notXUVar);

        // debug
//        cout << "Z2: ";
//        SymbolicSet Z2ss(*((*this->prodsX2s_[curAux])[curAbs]));
//        Z2ss.symbolicSet_ = Z2;
//        Z2ss.printInfo(1);

//        cout << "nZ2: ";
//        SymbolicSet nZ2ss(Z2ss);
//        nZ2ss.addGridPoints();
//        nZ2ss.symbolicSet_ &= nZ2;
//        nZ2ss.printInfo(1);

//        cout << "nC: ";
//        SymbolicSet nCss(*((*this->prodsCs_[curAux])[curAbs]));
//        nCss.addGridPoints();
//        nCss.symbolicSet_ &= nC;
//        nCss.printInfo(1);

//        cout << "C: ";
//        SymbolicSet Css(*((*this->prodsCs_[curAux])[curAbs]));
////        Css.addGridPoints();
//        Css.symbolicSet_ = C;
//        Css.printInfo(1);

        return C;
    }

    void innerFiner(SymbolicSet* Zc, SymbolicSet* Zf, SymbolicSet* XXc, BDD* notXVarf, SymbolicSet* Xf) {
        BDD Q = XXc->symbolicSet_ & Zc->symbolicSet_;
        Zf->symbolicSet_ = Q.ExistAbstract(*notXVarf) & Xf->symbolicSet_;
    }

    int innerCoarser(SymbolicSet* Zc, SymbolicSet* Zf, SymbolicSet* XXc, BDD* notXVarc, SymbolicSet* Xc, int numFiner) {
        // debug
//        XXc->printInfo(1);

        SymbolicSet Qcf(*XXc);
        Qcf.symbolicSet_ = XXc->symbolicSet_ & Zf->symbolicSet_;

        SymbolicSet Qc(*Zc);
        Qc.symbolicSet_ = Qcf.symbolicSet_.ExistAbstract(*notXVarc); // & S1
        Qc.symbolicSet_ &= !(Zc->symbolicSet_); /* don't check states that are already in Zc */

        int* QcMintermWhole;
        SymbolicSet Ccf(Qcf);
        BDD result = ddmgr_->bddZero();
        int* QcMinterm = new int[Qc.nvars_];

        for (Qc.begin(); !Qc.done(); Qc.next()) { // iterate over all coarse cells with any corresponding finer cells in Zf
            QcMintermWhole = (int*) Qc.currentMinterm();
            std::copy(QcMintermWhole + Qc.idBddVars_[0], QcMintermWhole + Qc.idBddVars_[0] + Qc.nvars_, QcMinterm);

            BDD coarseCell = Qc.mintermToBDD(QcMinterm) & Xc->symbolicSet_; // a particular coarse cell

            Ccf.symbolicSet_ = Qcf.symbolicSet_ & coarseCell; // corresponding finer cells to the coarse cell

            if ((Ccf.symbolicSet_.CountMinterm(Ccf.nvars_)) == (numFiner)) { // if there's a full set of finer cells
                result |= coarseCell;
            }
        }

        delete[] QcMinterm;

        if (result == ddmgr_->bddZero()) {
            return 0;
        }

        Zc->symbolicSet_ |= result;
        return 1;
    }

    void saveCZ(int curAux, int curAbs) {
        SymbolicSet* prodFinalC = new SymbolicSet(*((*this->prodsCs_[curAux])[curAbs]));
        SymbolicSet* prodFinalZ = new SymbolicSet(*((*this->prodsZs_[curAux])[curAbs]));
        prodFinalC->symbolicSet_ = ((*this->prodsCs_[curAux])[curAbs])->symbolicSet_;
        prodFinalZ->symbolicSet_ = ((*this->prodsZs_[curAux])[curAbs])->symbolicSet_;
        prodsFinalCs_[curAux]->push_back(prodFinalC);
        prodsFinalZs_[curAux]->push_back(prodFinalZ);
    }


    int reachOne(int startAbs, int minToGoCoarser, int minToBeValid, int earlyBreak, int verbose = 1) {
        startAbs_ = startAbs;
        minToGoCoarser_ = minToGoCoarser;
        minToBeValid_ = minToBeValid;
        earlyBreak_ = earlyBreak;
        verbose_ = verbose;

        if ( (startAbs_ < 0) || (startAbs_ >= *this->base_->numAbs_) ) {
            error("Error: startAbs out of range.\n");
        }
        if (this->auxs_.size() != 1) {
            error("Error: this debugging function only works with one auxiliary system.");
        }

        TicToc tt;
        tt.tic();

        int curAux = 0;
        int curAbs = startAbs_;
        int iter = 1;
        int justCoarsed = 0;
        int iterCurAbs = 1;
        int reached = 0;
        int stop = 0;

        // initialization
        for (int iAbs = 0; iAbs < *this->base_->numAbs_; iAbs++) {
            ((*this->prodsPreYs_[curAux])[iAbs])->addGridPoints();
            ((*this->prodsZs_[curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
            ((*this->prodsValidZs_[curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
            ((*this->prodsCs_[curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
            ((*this->prodsValidCs_[curAux])[iAbs])->symbolicSet_ = this->ddmgr_->bddZero();
        }

        while (1) {
            mu(curAux, &curAbs, &iter, &justCoarsed, &iterCurAbs, &reached, &stop);
            if (stop) {
                break;
            }
        }
        if (reached) {
            clog << "Won.\n";
            checkMakeDir("C");
            string prefix = "C/C";
            prefix += std::to_string(curAux+1);
            saveVec(*prodsFinalCs_[curAux], prefix);
            checkMakeDir("Z");
            prefix = "Z/Z";
            prefix += std::to_string(curAux+1);
            saveVec(*prodsFinalZs_[curAux], prefix);
            checkMakeDir("G");
            prefix = "G/G";
            prefix += std::to_string(curAux+1);
            saveVec(*prodsGs_[curAux], prefix);

            int fAbs = *this->base_->numAbs_ - 1;
            this->baseXs_[curAux]->writeToFile("plotting/X.bdd");
            this->baseOs_[fAbs]->writeToFile("plotting/O.bdd");
            ((*this->prodsIs_[curAux])[fAbs])->writeToFile("plotting/I.bdd");

            clog << "----------------------------------------reach: ";
            tt.toc();
            return 1;
        }
        else {
            clog << "Lost.\n";
            clog << "----------------------------------------reach: ";
            tt.toc();
            return 0;
        }
    }


};

}

#endif /* GBUCHI_HH_ */
