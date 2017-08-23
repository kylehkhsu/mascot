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

    int reachOne(int startAbs, int minToGoCoarser, int minToBeValid, int earlyBreak, int verbose = 1) {
        startAbs_ = startAbs;
        minToGoCoarser_ = minToGoCoarser;
        minToBeValid_ = minToBeValid;
        earlyBreak_ = earlyBreak;
        verbose_ = verbose;

        if ( (startAbs_ < 0) || (startAbs_ >= *this->base_->numAbs_) ) {
            error("Error: startAbs out of range.\n");
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


    void mu(int curAux, int* curAbs, int* iter, int* justCoarsed, int* iterCurAbs, int* reached, int* stop) {
        clog << "current abstraction: " << *curAbs << '\n';
        clog << "mu iteration: " << *iter << '\n';
        clog << "justCoarsed: " << *justCoarsed << '\n';
        clog << "iterCurAbs: " << *iterCurAbs << '\n';
        clog << "reached: " << *reached << '\n';
        clog << "controllers: " << this->prodsFinalCs_[curAux]->size() << '\n';

        BDD C = Cpre(((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_,
                ((*this->prodsTs_[curAux])[*curAbs])->symbolicSet_,
                ((*this->prodsTTs_[curAux])[*curAbs])->symbolicSet_,
                ((*this->prodsPermutesXtoX2_[curAux])[*curAbs]),
                *((*this->prodsNotXUVars_[curAux])[*curAbs]), curAux, *curAbs);
        C |= ((*this->prodsGs_[curAux])[*curAbs])->symbolicSet_;
        BDD N = C & (!((*this->prodsCs_[curAux])[*curAbs])->symbolicSet_.ExistAbstract(*((*this->prodsNotXVars_[curAux])[*curAbs])));
        ((*this->prodsCs_[curAux])[*curAbs])->symbolicSet_ |= N;
        ((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_ = C.ExistAbstract(*((*this->prodsNotXVars_[curAux])[*curAbs]));

        if ( (((*this->prodsZs_[curAux])[*curAbs])->symbolicSet_ & ((*this->prodsIs_[curAux])[*curAbs])->symbolicSet_) != this->ddmgr_->bddZero() ) {
            *reached = 1;
            if (earlyBreak_ == 1) {
                saveCZ(curAux, *curAbs);
                *stop = 1;
                clog << "\nTotal number of controllers: " << this->prodsFinalCs_[curAux]->size() << '\n';
                return;
            }
        }

        if (*iter != 1) {
            if (N == this->ddmgr_->bddZero()) {
                if (*curAbs == *this->base_->numAbs_ - 1) {
                    if (*reached == 1) {
                        saveCZ(curAux, *curAbs);
                    }
                    *stop = 1;
                    clog << "\nTotal number of controllers: " << this->prodsFinalCs_[curAux]->size() << '\n';
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

    BDD Cpre (BDD Z, BDD T, BDD TT, int* permuteXtoX2, BDD notXUVar, int curAux, int curAbs) {
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
        cout << "0\n";
        // debug
//        XXc->printInfo(1);


        SymbolicSet Qcf(*XXc);
        Qcf.symbolicSet_ = XXc->symbolicSet_ & Zf->symbolicSet_;

        cout << "1\n";
        SymbolicSet Qc(*Zc);
        Qc.symbolicSet_ = Qcf.symbolicSet_.ExistAbstract(*notXVarc); // & S1
        Qc.symbolicSet_ &= !(Zc->symbolicSet_); /* don't check states that are already in Zc */

        int* QcMintermWhole;
        SymbolicSet Ccf(Qcf);
        BDD result = ddmgr_->bddZero();
        int* QcMinterm = new int[Qc.nvars_];

        cout << "2\n";
        for (Qc.begin(); !Qc.done(); Qc.next()) { // iterate over all coarse cells with any corresponding finer cells in Zf
            QcMintermWhole = (int*) Qc.currentMinterm();
            std::copy(QcMintermWhole + Qc.idBddVars_[0], QcMintermWhole + Qc.idBddVars_[0] + Qc.nvars_, QcMinterm);

            BDD coarseCell = Qc.mintermToBDD(QcMinterm) & Xc->symbolicSet_; // a particular coarse cell

            Ccf.symbolicSet_ = Qcf.symbolicSet_ & coarseCell; // corresponding finer cells to the coarse cell

            if ((Ccf.symbolicSet_.CountMinterm(Ccf.nvars_)) == (numFiner)) { // if there's a full set of finer cells
                result |= coarseCell;
            }
        }
        cout << "3\n";
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


};

}

#endif /* GBUCHI_HH_ */
