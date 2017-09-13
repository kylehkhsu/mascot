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
//                cout << "Z:\n";
//                this->Zs_[curAbs]->printInfo(1);
                break;
            }
            iter += 1;
        }

        infZs_[curAbs]->symbolicSet_ |= Zs_[curAbs]->symbolicSet_;
        if (curAbs != *this->system_->numAbs_ - 1) {
            this->finer(this->infZs_[curAbs], this->infZs_[curAbs+1], curAbs); // obtain projection of converged Z onto next, finer abstraction
//            cout << "infZ:\n";
//            this->infZs_[curAbs+1]->printInfo(1);

            this->Ts_[curAbs+1]->symbolicSet_ &= !(infZs_[curAbs+1]->symbolicSet_); // don't consider pre states that already have a controller in a coarser abstraction
            *this->TTs_[curAbs+1] &= !(infZs_[curAbs+1]->symbolicSet_); // same as above
        }
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
        this->infZs_[*this->system_->numAbs_-1]->printInfo(1);
        this->infZs_[*this->system_->numAbs_-1]->writeToFile("plotting/D.bdd");
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
