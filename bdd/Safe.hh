/*! \file Safe.hh
 *  Contains the Safe class. */

#ifndef SAFE_HH_
#define SAFE_HH_

#include "Adaptive.hh"

using namespace helper;

namespace scots {

/*! \class Safe
 *  \brief A class (derived from base Adaptive) that does adaptive multiscale abstraction-based synthesis for a safe specification.
 */
class Safe: public virtual Adaptive {
public:
    vector<SymbolicSet*> Ss_; /*!< Instance of *Xs_[i] containing safe states. */
    vector<SymbolicSet*> infZs_; /*!< Instance of *Xs_[i] containing projection of convergence of previous maximal fixed points. */
    vector<SymbolicSet*> validZs_; /*!< Contains winning states that act as savepoints. */


    /*! Constructor for a Safe object. */
    Safe(char* logFile)
        : Adaptive(logFile)
    {
    }

    ~Safe() {
        deleteVec(Ss_);
        deleteVec(infZs_);
        deleteVec(validZs_);
    }

    void muSafe(int minToGoCoarser, int minToBeValid, int verbose, int* curAbs, int* iter, int* justCoarsed, int* iterCurAbs, int* stop) {
        clog << "current abstraction: " << *curAbs << '\n';
        clog << "mu iteration: " << *iter << '\n';
        clog << "justCoarsed: " << *justCoarsed << '\n';
        clog << "iterCurAbs: " << *iterCurAbs << '\n';

        cout << "current abstraction: " << *curAbs << '\n';
        cout << "mu iteration: " << *iter << '\n';
        cout << "justCoarsed: " << *justCoarsed << '\n';
        cout << "iterCurAbs: " << *iterCurAbs << '\n';

//        // get pre(Z)
//        BDD nPreZ = this->nPreC(this->Zs_[*curAbs]->symbolicSet_, this->Ts_[*curAbs]->symbolicSet_, *this->TTs_[*curAbs], *curAbs);
//        // disjunct with safe set (minimal fixed point)
//        nPreZ |= !Ss_[*curAbs]->symbolicSet_;

        if (*iter == 50) {
            *stop = 1;
        }


        cout << "Z:\n";
        this->Zs_[*curAbs]->printInfo(1);

        BDD nZ = !this->Zs_[*curAbs]->symbolicSet_ & this->Xs_[*curAbs]->symbolicSet_;

        SymbolicSet nZSS(*this->Zs_[*curAbs]);
        nZSS.symbolicSet_ = nZ;
//        cout << "nZ:\n";
//        nZSS.printInfo(1);

        BDD preC = this->preC(nZSS.symbolicSet_, this->Ts_[*curAbs]->symbolicSet_, *this->TTs_[*curAbs], *curAbs) & this->Xs_[*curAbs]->symbolicSet_ & this->U_->symbolicSet_;

//        SymbolicSet preCSS(*this->Cs_[*curAbs]);
//        preCSS.symbolicSet_ = preC;
//        cout << "preC:\n";
//        preCSS.printInfo(1);

        BDD preZ = preC.ExistAbstract(*notXvars_[*curAbs]);

        SymbolicSet preZSS(*this->Zs_[*curAbs]);
//        preZSS.symbolicSet_ = preZ;
//        cout << "preZ:\n";
//        preZSS.printInfo(1);

        preZ = ((!preZ) | (!Ss_[*curAbs]->symbolicSet_)) & this->Xs_[*curAbs]->symbolicSet_;

        preZSS.symbolicSet_ = preZ;
        cout << "preZ:\n";
        preZSS.printInfo(1);

        if (preZ == Zs_[*curAbs]->symbolicSet_) { // if iteration failed to yield new (x,u)
            if (*curAbs == *this->system_->numAbs_ - 1) { // if we're in the finest abstraction
                *stop = 1;
            }
            else {
                if (verbose) {
                    clog << "No new winning states; going finer.\n";
                }
                if (*justCoarsed == 1) { // current winning set has not been declared valid
                    if (verbose) {
                        clog << "Current winning set has not been declared valid.\n";
                        clog << "Resetting Zs_[" << *curAbs << "] from valids.\n";
                    }

                    // reset the current winning states and controller
                    this->Zs_[*curAbs]->symbolicSet_ = validZs_[*curAbs]->symbolicSet_;
                }
                else {
                    if (verbose) {
                        clog << "Current winning set has been declared valid.\n";
                        clog << "Saving Z of abstraction " << *curAbs << " into valids.\n";
                        clog << "Finding inner approximation of Zs_[" << *curAbs << "] and copying into validZs_[" << *curAbs+1 << "].\n";
                    }
                    validZs_[*curAbs]->symbolicSet_ = this->Zs_[*curAbs]->symbolicSet_;

                    this->innerFinerAligned(this->Zs_[*curAbs], this->Zs_[*curAbs+1], *curAbs);
                    validZs_[*curAbs+1]->symbolicSet_ = this->Zs_[*curAbs+1]->symbolicSet_;
                }
                *curAbs += 1;
                *iterCurAbs = 1;
            }
        }
        else { // if there were new (x,u)
            Zs_[*curAbs]->symbolicSet_ = preZ;
            if ((*justCoarsed == 1) && (*iterCurAbs >= minToBeValid)) {
                if (verbose) {
                    clog << "Current winning set now valid.\n";
                }
                *justCoarsed = 0;
            }
            if (*curAbs == 0) {
                if (verbose) {
                    clog << "More new states, already at coarsest gridding, continuing.\n";
                }
                *iterCurAbs += 1;
            }
            else {
                if (*iterCurAbs >= minToGoCoarser) {
                    if (verbose) {
                        clog << "More new states, minToGoCoarser achieved; try going coarser.\n";
                    }
                    int more = this->innerCoarserAligned(this->Zs_[*curAbs-1], this->Zs_[*curAbs], *curAbs-1);

                    if (more == 0) {
                        if (verbose) {
                            clog << "Projecting to coarser gives no more states in coarser; not changing abstraction.\n";
                        }
                        *iterCurAbs += 1;
                    }
                    else {
                        if (verbose) {
                            clog << "Projecting to coarser gives more states in coarser.\n";
                            clog << "Saving Z of abstraction " << *curAbs << " into validZs.\n";
                            clog << "Saving Z of abstraction " << *curAbs-1 << " into validZs_.\n";
                        }
                        validZs_[*curAbs]->symbolicSet_ = this->Zs_[*curAbs]->symbolicSet_;
                        *justCoarsed = 1;
                        validZs_[*curAbs-1]->symbolicSet_ = this->Zs_[*curAbs-1]->symbolicSet_;
                        *curAbs -= 1;
                        *iterCurAbs = 1;
                    }
                }
                else {
                    *iterCurAbs += 1;
                }
            }
        }
        clog << "\n";
        *iter += 1;
    }

    /*! Adaptive maximal fixed point. */
    void nu(int curAbs) {
        int iter = 1;
        while (1) {
            clog << ".";
            // get pre of current abtraction's Z disjuncted with projection of converged Z from previous abstraction
            this->Cs_[curAbs]->symbolicSet_ = this->preC(this->Zs_[curAbs]->symbolicSet_ | infZs_[curAbs]->symbolicSet_, this->Ts_[curAbs]->symbolicSet_, *this->TTs_[curAbs], curAbs);
//            this->Cs_[curAbs]->symbolicSet_ = this->preC(this->Zs_[curAbs]->symbolicSet_, this->Ts_[curAbs]->symbolicSet_, *this->TTs_[curAbs], curAbs);
            // conjunct with safe set (maximal fixed point)
//            this->Cs_[curAbs]->symbolicSet_ &= (Ss_[curAbs]->symbolicSet_ & !infZs_[curAbs]->symbolicSet_);
            this->Cs_[curAbs]->symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
            // project onto Xs_[curAbs]
            BDD Z = this->Cs_[curAbs]->symbolicSet_.ExistAbstract(*this->notXvars_[curAbs]);

            if (Z != this->Zs_[curAbs]->symbolicSet_) { // not converged
                this->Zs_[curAbs]->symbolicSet_ = Z; // update and continue
            }
            else { // converged
                clog << iter << '\n';
                this->Zs_[curAbs]->printInfo(1);
                break;
            }
            iter += 1;
        }
        if (curAbs != *this->system_->numAbs_ - 1) {
            this->innerFinerAligned(this->Zs_[curAbs], this->infZs_[curAbs+1], curAbs); // obtain projection of converged Z onto next, finer abstraction

            this->Ts_[curAbs+1]->symbolicSet_ &= !(infZs_[curAbs+1]->symbolicSet_); // don't consider pre states that already have a controller in a coarser abstraction
            *this->TTs_[curAbs+1] &= !(infZs_[curAbs+1]->symbolicSet_); // same as above
        }
        this->Xs_[curAbs]->symbolicSet_ = this->Zs_[curAbs]->symbolicSet_ | infZs_[curAbs]->symbolicSet_; // for verification purposes
    }

    void safe2(int startAbs, int minToGoCoarser, int minToBeValid, int verbose = 1) {
        TicToc tt;
        tt.tic();

        // removing obstacles from transition relation
        this->removeFromTs(&(this->Os_));

        clog << "Os_ removed from Ts_, TTs_.\n";

        int curAbs = startAbs;
        int iter = 1;
        int justCoarsed = 0;
        int iterCurAbs = 1;
        int stop = 0;

        while (1) {
            muSafe(minToGoCoarser, minToBeValid, verbose, &curAbs, &iter, &justCoarsed, &iterCurAbs, &stop);
            if (stop) {
                break;
            }
        }

        infZs_[0]->symbolicSet_ = !this->Zs_[0]->symbolicSet_;
        for (int i = 0; i < *this->system_->numAbs_-1; i++) {
            this->innerFinerAligned(infZs_[i], infZs_[i+1], i);
            infZs_[i+1]->symbolicSet_ &= !this->Zs_[i+1]->symbolicSet_;
        }

        int fAb = *this->system_->numAbs_ - 1;

//        cout << "infZ:\n";
//        infZs_[fAb]->symbolicSet_ &= Xs_[fAb]->symbolicSet_;
//        infZs_[fAb]->printInfo(1);

//        cout << "Z:\n";
//        Zs_[fAb]->symbolicSet_ &= Xs_[fAb]->symbolicSet_;
//        Zs_[fAb]->printInfo(1);

        this->Zs_[fAb]->symbolicSet_ = infZs_[fAb]->symbolicSet_ & Xs_[fAb]->symbolicSet_;
        cout << "Z:\n";
        this->Zs_[fAb]->printInfo(1);

        this->Cs_[fAb]->symbolicSet_ = this->preC(this->Zs_[fAb]->symbolicSet_, this->Ts_[fAb]->symbolicSet_, *this->TTs_[fAb], fAb);

        cout << "C finest:\n";
        Cs_[fAb]->printInfo(1);

        SymbolicSet lol(*Zs_[fAb]);
        lol.symbolicSet_ = Cs_[fAb]->symbolicSet_.ExistAbstract(*notXvars_[fAb]);
        cout << "Z again:\n";
        lol.printInfo(1);

        Cs_[fAb]->writeToFile("C/C.bdd");


    }

    /*! Writes, should they exist, controllers of specified abstractions that together satisfy a safety specification. */
    void safe() {

        TicToc tt;
        tt.tic();

        // removing obstacles from transition relation
        this->removeFromTs(&(this->Os_));

        for (int curAbs = 0; curAbs < *this->system_->numAbs_; curAbs++) {
            nu(curAbs);
        }

        clog << "----------------------------------------safe: ";
        tt.toc();

        checkMakeDir("Z");
        saveVec(this->Zs_, "Z/Z");
        checkMakeDir("C");
        saveVec(this->Cs_, "C/C");

        cout << "Domain of all controllers:\n";
        this->Xs_[*this->system_->numAbs_-1]->printInfo(1);
        this->Xs_[*this->system_->numAbs_-1]->writeToFile("plotting/D.bdd");
    }

    /*! Initializes objects specific to the following specifications: safe.
        \param[in]	addS	Function pointer specifying the points that should be added to the potential safe set.
    */
    template<class S_type>
    void initializeSafe(S_type addS) {

        TicToc tt;
        tt.tic();
        this->initializeSpec(&Ss_, addS);
        clog << "Ss_ initialized\n";
        clog << "------------------------------------------------initializeSpec on Ss_: ";
        tt.toc();

        // maximal fixed point starts with whole set
        for (int i = 0; i < *this->system_->numAbs_; i++) {
            this->Zs_[i]->symbolicSet_ = this->ddmgr_->bddOne();
            ;
        }
        clog << "Zs_ modified to BDD-1.\n";

        for (int i = 0; i < *this->system_->numAbs_; i++) {
            SymbolicSet* infZ = new SymbolicSet(*this->Xs_[i]);
            infZs_.push_back(infZ);
        }
        clog << "infZs_ initialized to empty.\n";

        for (int i = 0; i < *this->system_->numAbs_; i++) {
            SymbolicSet* validZ = new SymbolicSet(*this->Xs_[i]);
            validZs_.push_back(validZ);
        }
        clog << "validZs_ initialized with empty domain.\n";

        saveVerify();
    }

    /*! Saves and prints to log file some information related to the safety specification. */
    void saveVerify() {
        printVec(Ss_, "S");
        checkMakeDir("S");
        saveVec(Ss_,"S/S");
    }

    //    /*! Debugging fucntion. Implementation of basic SCOTS safety (single abstraction) in the Adaptive framework. */
    //    void safeBasicDebug(int verbose = 0) {
    //        if (numAbs_ != 1) {
    //            error("Error: For comparison with SCOTS, numAbs needs to be 1.\n");
    //        }

    //        int curAbs = 0;

    //        TicToc tt;
    //        tt.tic();

    //        SymbolicSet C(*Xs_[curAbs], *U_);

    //        int iter = 1;

    //        while (1) {
    //            // get pre(Z)
    //            C.symbolicSet_ = preC(Zs_[curAbs]->symbolicSet_, curAbs);
    //            // conjunct with safe set (maximal fixed point)
    //            C.symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
    //            // project onto Xs_[curAbs]
    //            BDD Z = C.symbolicSet_.ExistAbstract(*notXvars_[curAbs]);

    //            clog << "iter: " << iter << '\n';
    //            iter += 1;

    //            if (Z != Zs_[curAbs]->symbolicSet_) {
    //                Zs_[curAbs]->symbolicSet_ = Z;
    //            }
    //            else {
    //                break;
    //            }
    //        }

    //        tt.toc();

    //        C.printInfo(1);
    //        checkMakeDir("debug");
    //        C.writeToFile("debug/C.bdd");

    //        SymbolicSet X(*Xs_[curAbs]);
    //        X.symbolicSet_ = C.symbolicSet_.ExistAbstract(U_->getCube());

    //        clog << "Domain of controller:\n";
    //        X.printInfo(1);
    //    }

};
}

#endif /* SAFE_HH_ */
