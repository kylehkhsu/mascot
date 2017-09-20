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

    int itersToNextAbs_;

    /*! Constructor for a Safe object. */
    Safe(char* logFile)
        : Adaptive(logFile)
    {
    }

    ~Safe() {
        deleteVec(Ss_);
        deleteVec(infZs_);
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
            this->Cs_[curAbs]->symbolicSet_ &= (Ss_[curAbs]->symbolicSet_ & !infZs_[curAbs]->symbolicSet_);
//            this->Cs_[curAbs]->symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
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
        infZs_[curAbs]->printInfo(1);
        if (curAbs != *this->system_->numAbs_ - 1) {
            this->finer(infZs_[curAbs], infZs_[curAbs+1], curAbs); // obtain projection of converged Z onto next, finer abstraction

//            this->Ts_[curAbs+1]->symbolicSet_ &= !(infZs_[curAbs+1]->symbolicSet_); // don't consider pre states that already have a controller in a coarser abstraction
//            *this->TTs_[curAbs+1] &= !(infZs_[curAbs+1]->symbolicSet_); // same as above
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
        infZs_[*this->system_->numAbs_-1]->printInfo(1);
        infZs_[*this->system_->numAbs_-1]->writeToFile("plotting/D.bdd");
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
        }
        clog << "Zs_ modified to BDD-1.\n";

        for (int i = 0; i < *this->system_->numAbs_; i++) {
            SymbolicSet* infZ = new SymbolicSet(*this->Xs_[i]);
            infZs_.push_back(infZ);
        }
        clog << "infZs_ initialized to empty.\n";

        saveVerify();
    }

    /*! Saves and prints to log file some information related to the safety specification. */
    void saveVerify() {
        printVec(Ss_, "S");
        checkMakeDir("S");
        saveVec(Ss_,"S/S");
    }

    void nu2(int* curAbs, int* iterTotal, int* stop) {
        int iterCurAbs = 1;
        clog << "current abstraction: " << *curAbs << '\n';
        while (iterCurAbs <= itersToNextAbs_) {
            *iterTotal += 1;
            clog << ".";
            // get pre of current abtraction's Z disjuncted with projection of Z from previous abstraction
            this->Cs_[*curAbs]->symbolicSet_ = this->preC(this->Zs_[*curAbs]->symbolicSet_ | infZs_[*curAbs]->symbolicSet_, this->Ts_[*curAbs]->symbolicSet_, *this->TTs_[*curAbs], *curAbs);
            // conjunct with safe set (maximal fixed point)
            this->Cs_[*curAbs]->symbolicSet_ &= (Ss_[*curAbs]->symbolicSet_ & !infZs_[*curAbs]->symbolicSet_);
//            this->Cs_[*curAbs]->symbolicSet_ &= Ss_[*curAbs]->symbolicSet_;
            // project onto Xs_[curAbs]
            BDD Z = this->Cs_[*curAbs]->symbolicSet_.ExistAbstract(*this->notXvars_[*curAbs]);
            if (Z != this->Zs_[*curAbs]->symbolicSet_) { // not converged
                this->Zs_[*curAbs]->symbolicSet_ = Z; // update and continue
                iterCurAbs += 1;
            }
            else { // this abstraction has converged
                if (*curAbs == *this->system_->numAbs_ - 1) { // this is the finest abstraction
                    *stop = 1;
                }
                break;
            }
        }
        clog << iterCurAbs << " ";
        clog << "total iterations: " << *iterTotal << '\n';
        infZs_[*curAbs]->symbolicSet_ |= this->Zs_[*curAbs]->symbolicSet_;
        infZs_[*curAbs]->printInfo(1);
        if (*this->system_->numAbs_ == 1) {
            return;
        }
        else {
            if (*curAbs == *this->system_->numAbs_ - 1) { // finished a pass
                if (*stop == 1) {
                    return;
                }
                else {
                    for (int iAbs = *curAbs; iAbs > 0; iAbs--) {
                        this->coarser(infZs_[iAbs - 1], infZs_[iAbs], iAbs - 1, 0);
                    }
                    this->Ts_[0]->symbolicSet_ &= !(infZs_[0]->symbolicSet_);
                    *this->TTs_[0] &= !(infZs_[0]->symbolicSet_);

                    for (int i = 0; i < *this->system_->numAbs_; i++) {
                        this->Zs_[i]->symbolicSet_ = infZs_[i]->symbolicSet_;
                    }
                    clog << "Zs updated to infZs.\n";

                    *curAbs = 0;
                }
            }
            else {
                this->finer(infZs_[*curAbs], infZs_[*curAbs+1],  *curAbs);
//                this->Ts_[*curAbs+1]->symbolicSet_ &= !(infZs_[*curAbs+1]->symbolicSet_); // don't consider pre states that already have a controller in a coarser abstraction
//                *this->TTs_[*curAbs+1] &= !(infZs_[*curAbs+1]->symbolicSet_); // same as above
                *curAbs += 1;
            }
        }
    }

    void safe2(int itersToNextAbs) {

        itersToNextAbs_ = itersToNextAbs;
        TicToc tt;
        tt.tic();

        // removing obstacles from transition relation
        this->removeFromTs(&(this->Os_));
        clog << "Os_ removed from Ts_, TTs_.\n";


        int curAbs = 0;
        int iterTotal = 0;
        int stop = 0;
        while (1) {
            nu2(&curAbs, &iterTotal, &stop);
            if (stop == 1) {
                break;
            }
        }

        clog << "----------------------------------------safe: ";
        tt.toc();

        checkMakeDir("Z");
        saveVec(this->Zs_, "Z/Z");
        checkMakeDir("C");
        saveVec(this->Cs_, "C/C");

        cout << "Domain of all controllers:\n";
        infZs_[*this->system_->numAbs_-1]->printInfo(1);
        infZs_[*this->system_->numAbs_-1]->writeToFile("plotting/D.bdd");

    }
};
}

#endif /* SAFE_HH_ */
