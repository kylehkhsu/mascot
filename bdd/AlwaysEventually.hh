/*! \file AlwaysEventually.hh
 *  Contains the AlwaysEventually class. */

#ifndef ALWAYSEVENTUALLY_HH_
#define ALWAYSEVENTUALLY_HH_

#include "Reach.hh"

using namespace helper;

namespace scots {

/*! \class AlwaysEventually
 *  \brief A class (derived from base Adaptive) that does adaptive multiscale abstraction-based synthesis for an always-eventually specification.
 */
template<class X_type, class U_type>
class AlwaysEventually: public virtual Reach<X_type, U_type> {
public:

    /*! Constructor for an AlwaysEventually object. */
    AlwaysEventually(int dimX, double* lbX, double* ubX, double* etaX, double tau,
                     int dimU, double* lbU, double* ubU, double* etaU,
                     double* etaRatio, double tauRatio, int nint,
                     int numAbs, int readXX, int readAbs, char* logFile)
                     : Reach<X_type, U_type>(dimX, lbX, ubX, etaX, tau,
                                             dimU, lbU, ubU, etaU,
                                             etaRatio, tauRatio, nint,
                                             numAbs, readXX, readAbs, logFile),
                       Adaptive<X_type, U_type>(dimX, lbX, ubX, etaX, tau,
                                                dimU, lbU, ubU, etaU,
                                                etaRatio, tauRatio, nint,
                                                numAbs, readXX, readAbs, logFile)
    {
    }

    /*! Deconstructor for an AlwaysEventually object. */
    ~AlwaysEventually() {

    }

    /*!	Writes, should they exist, a sequence of controller and controller domain BDDs to directories 'C' and 'Z' respectively that satisfy the Buchi box-diamond (aka always-eventually) specification.
     *  Note: of the n resultant controllers in 'C', the n-1th, .., 1st controllers are to reach, and the nth controller is for when the system state is a goal state.
        \param[in]  startAbs            0-index of the abstraction to start with.
        \param[in]	minToGoCoarser		Minimum number of growing fixed point iterations needed before an attempt to go to a coarser abstraction.
        \param[in]	minToBeValid		Minimum number of growing fixed point iterations needed before a controller is declared valid.
        \param[in]	verbose				If 1, prints additional information during synthesis to the log file.
        \return     1 if controller(s) satisfying specification is/are synthesized; 0 otherwise.
     */
    int alwaysEventually(int startAbs, int minToGoCoarser, int minToBeValid, int verbose = 1) {
        if (this->stage_ != 3) {
            error("Error: alwaysEventually called out of order.\n");
        }
        if ( (startAbs < 0) || (startAbs >= this->numAbs_) ) {
            error("Error: startAbs out of range.\n");
        }
        int earlyBreak = 0;

        TicToc tt;
        tt.tic();

        // removing obstacles from transition relation
        this->removeFromTs(&(this->Os_));

        clog << "Os_ removed from Ts_, TTs_.\n";


        int curAbs = startAbs;

        int nuIter = 1;

        while (1) {
            clog << "------------------------------nu iteration: " << nuIter << "------------------------------\n";

            int muIter = 1;
            int justCoarsed = 0;
            int iterCurAbs = 1;
            int reached = 0;
            int muStop = 0;
            while (1) { // eventually
                this->mu(minToGoCoarser, minToBeValid, earlyBreak, verbose, &curAbs, &muIter, &justCoarsed, &iterCurAbs, &reached, &muStop);
                if (muStop) {
                    break;
                }
            }

            clog << "\nAlwaysEventually: iteration " << nuIter << "\n";

            // always
            // project to the finest abstraction
            for (int i = curAbs; i < this->numAbs_ - 1; i++) {
                this->innerFinerAligned(this->Zs_[i], this->Zs_[i+1], i);
            }
            this->Zs_[this->numAbs_ - 1]->symbolicSet_ &= this->Xs_[this->numAbs_ - 1]->symbolicSet_;
            curAbs = this->numAbs_ - 1;


            BDD preQ = this->preC(this->Zs_[curAbs]->symbolicSet_, this->Ts_[curAbs]->symbolicSet_, *this->TTs_[curAbs], curAbs); // preC(X2)
            BDD thing = preQ & this->Gs_[curAbs]->symbolicSet_; // R \cap preC(X2); converges to controller for goal states
            BDD Q = thing.ExistAbstract(*this->notXvars_[curAbs]); // converges to goal states that are valid for the AE spec (i.e. can continue infinitely)

            if ((this->Gs_[curAbs]->symbolicSet_ <= Q) && (Q <= this->Gs_[curAbs]->symbolicSet_)) {
                clog << "Won.\n";

                SymbolicSet* lastZ = this->finalZs_.back();
                int abs = -1;
                for (int i = 0; i < this->numAbs_; i++) {
                    int sameEta = 1;
                    for (int j = 0; j < this->dimX_; j++) {
                        if (this->etaX_[i][j] != lastZ->eta_[j]) {
                            sameEta = 0;
                        }
                    }
                    int sameTau = this->tau_[i][0] == lastZ->tau_;

                    if (sameEta && sameTau) {
                        abs = i;
                        break;
                    }
                }

                clog << "Abstraction of E: " << abs << '\n';

                SymbolicSet E(*this->Xs_[abs]);
                E.symbolicSet_ = lastZ->symbolicSet_ & !(this->Gs_[abs]->symbolicSet_);
                E.writeToFile("G/E.bdd");

                SymbolicSet* goalC = new SymbolicSet(*this->Cs_[curAbs]);
                goalC->symbolicSet_ = thing;
                this->finalCs_.push_back(goalC);

                SymbolicSet* goalZ = new SymbolicSet(*this->Zs_[curAbs]);
                goalZ->symbolicSet_ = Q;
                this->finalZs_.push_back(goalZ);

                clog << this->finalCs_.size()<< '\n';

                checkMakeDir("C");
                saveVec(this->finalCs_, "C/C");
                checkMakeDir("Z");
                saveVec(this->finalZs_, "Z/Z");
                clog << "----------------------------------------alwaysEventually: ";
                tt.toc();

                return 1;
            }
            else {
                clog << "Continuing.\n";


                clog << "Updating Gs.\n";
                this->Gs_[curAbs]->symbolicSet_ = Q;
                for (int i = curAbs; i > 0; i--) {
                    this->innerCoarserAligned(this->Gs_[i-1], this->Gs_[i], i-1);
                }

                clog << "Resetting Zs, Cs.\n";
                for (int i = 0; i < this->numAbs_; i++) {
                    this->Zs_[i]->symbolicSet_ = this->ddmgr_->bddZero();
                    this->Cs_[i]->symbolicSet_ = this->ddmgr_->bddZero();
                }

                size_t j = this->finalCs_.size();
                for (size_t i = 0; i < j; i++) {
                    delete this->finalCs_.back();
                    this->finalCs_.pop_back();
                    delete this->finalZs_.back();
                    this->finalZs_.pop_back();
                }

                clog << "finalCs_.size(): " << this->finalCs_.size() << '\n';

            }
            clog << '\n';
            nuIter += 1;
        }
    }
};

}

#endif /* ALWAYSEVENTUALLY_HH_ */
