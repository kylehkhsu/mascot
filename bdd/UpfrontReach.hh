/*! \file Reach.hh
 *  Contains the Reach class. */

#ifndef UPFRONTREACH_HH_
#define UPFRONTREACH_HH_

#include "Adaptive.hh"

using namespace helper;

namespace scots {

enum ReachResult {CONVERGEDVALID, CONVERGEDINVALID, NOTCONVERGED};

/*! \class Reach
 *  \brief A class (derived from base Adaptive) that does adaptive multiscale abstraction-based synthesis for a reach specification.
 */
class UpfrontReach: public virtual Adaptive {
public:

    vector<SymbolicSet*> Gs_; /*!< Instance of *Xs_[i] containing goal states. */
    vector<SymbolicSet*> Is_; /*!< Instance of *Xs_[i] containing initial states. */
    vector<SymbolicSet*> validZs_; /*!< Contains winning states that act as savepoints. */
    vector<SymbolicSet*> validCs_; /*!< Controllers that act as savepoints. */
    vector<SymbolicSet*> finalCs_; /*!< Sequence of controllers that satisfy the specification. */
    vector<SymbolicSet*> finalZs_; /*!< Sequence of domains of finalCs_. */

    int m_;
    int earlyBreak_;
    int verbose_;


    /*! Constructor for a Reach object. */
    UpfrontReach(char* logFile)
        : Adaptive(logFile)
    {
    }

    /*! Deconstructor for a Reach object. */
    ~UpfrontReach() {
        deleteVec(Gs_);
        deleteVec(Is_);
        deleteVec(validZs_);
        deleteVec(validCs_);
        deleteVec(finalCs_);
        deleteVec(finalZs_);
    }

    /*! Calculates the enforceable predecessor of the given set with respect to the transition relation at the specified level of abstraction.
     *  \param[in]  Z           The winning set.
     *  \param[in]  curAbs      0-index of the current abstraction.
     *  \return     BDD containing {(x,u)} for which all post states are in Z.
     */
    BDD cPre(BDD Z, int curAbs) {
        // swap to X2s_[curAbs]
        BDD Z2 = Z.Permute(this->permutesXtoX2_[curAbs]);
        // posts outside of winning set
        BDD nZ2 = !Z2;
        // {(x,u)} with some posts outside of winning set
        BDD Fbdd = this->Ts_[curAbs]->symbolicSet_.AndAbstract(nZ2, *this->notXUvars_[curAbs]);
        // {(x,u)} with no posts outside of winning set
        BDD nF = !Fbdd;
        // get rid of junk
        BDD preF = this->TTs_[curAbs]->AndAbstract(nF, *this->notXUvars_[curAbs]);
        return preF;
    }

    ReachResult reach(int ab, int m = -1, int print = 0) {
        int i = 1;
        clog << "abstraction: " << ab << '\n';
        if (print)
            cout << "abstraction: " << ab << '\n';
        while (1) {
            clog << "iteration: " << i << '\n';
            if (print)
                cout << "iteration: " << i << '\n';

//            if (ab != 0) {
//                Zs_[ab]->printInfo(1);
//            }
            BDD cPreZ = cPre(Zs_[ab]->symbolicSet_, ab);
            BDD C = cPreZ | Gs_[ab]->symbolicSet_;
            BDD N = C & (!(Cs_[ab]->symbolicSet_.ExistAbstract(*this->cubeU_)));
            Cs_[ab]->symbolicSet_ |= N;
            Zs_[ab]->symbolicSet_ = C.ExistAbstract(*notXvars_[ab]) | validZs_[ab]->symbolicSet_; // the disjunction part is new

            if (N == ddmgr_->bddZero() && i != 1) {
                if (i >= m_) {
                    return CONVERGEDVALID;
                }
                else{
                    return CONVERGEDINVALID;
                }
            }
            if (i == m) {
                return NOTCONVERGED;
            }
            i += 1;
        }
    }

    void finer(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        Zf->addGridPoints();
        Zf->symbolicSet_ &= Zc->symbolicSet_.Permute(permutesFiner_[c]);
    }

    int coarserInner(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        BDD Zcand = Zf->symbolicSet_.UnivAbstract(*cubesCoarser_[c]);
        Zcand = Zcand.Permute(permutesCoarser_[c]);

        if (Zcand <= Zc->symbolicSet_) {
            return 0;
        }
        else {
            Zc->symbolicSet_ = Zcand;
            return 1;
        }
    }

    void coarserOuter(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        Zc->symbolicSet_ = (Zf->symbolicSet_.ExistAbstract(*cubesCoarser_[c])).Permute(permutesCoarser_[c]);
    }

    void reachRecurse(int ab) {
        clog << '\n';
        clog << "current abstraction: " << ab << '\n';
        clog << "controllers: " << finalCs_.size() << '\n';
        ReachResult result;
        if (ab == 0) {
            result = reach(ab);
        }
        else {
            result = reach(ab, m_);
        }
        if (result == CONVERGEDVALID) {
            clog << "result: converged valid\n";
        }
        else if (result == CONVERGEDINVALID) {
            clog << "result: converged invalid\n";
        }
        else {
            clog << "result: not converged\n";
        }

        if (result != CONVERGEDINVALID) {
            saveCZ(ab);
            validZs_[ab]->symbolicSet_ = Zs_[ab]->symbolicSet_;
            validCs_[ab]->symbolicSet_ = Cs_[ab]->symbolicSet_;
            clog << "saved as snapshot, saved to valids\n";
        }
        else { // result is CONVERGEDINVALID
            Zs_[ab]->symbolicSet_ = validZs_[ab]->symbolicSet_; // reset this layer's progress
            Cs_[ab]->symbolicSet_ = validCs_[ab]->symbolicSet_;
            clog << "reset to valids\n";
        }

        if (result != NOTCONVERGED) {
            if (ab == *system_->numAbs_ - 1) {
                return;
            }
            else { // go finer
                clog << "Going finer\n";
                int nextAb = ab + 1;
                finer(Zs_[ab], Zs_[nextAb], ab);
                Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
                validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
                reachRecurse(nextAb);
                return;
            }
        }
        else { // not converged, go coarser
            clog << "Going coarser\n";
            int nextAb = ab - 1;
            coarserInner(Zs_[nextAb], Zs_[ab], nextAb);
            Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
            validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
            reachRecurse(nextAb);
            return;
        }
    }

    void upfrontReach(int m, int earlyBreak = 0, int verbose = 1) {
        m_ = m;
        earlyBreak_ = earlyBreak;
        verbose_ = verbose;

        clog << "Os_ removed from Ts_, TTs_.\n";

        TicToc timer;
        timer.tic();
        int ab = 0;
        reachRecurse(ab);

        clog << "controllers: " << finalCs_.size() << '\n';

        checkMakeDir("C");
        saveVec(finalCs_, "C/C");
        checkMakeDir("Z");
        saveVec(finalZs_, "Z/Z");
        clog << "----------------------------------------Upfront Reach: " << timer.toc() << " seconds.\n";
    }

    /*! Initializes objects specific to the following specifications: always-eventually, reach-while-avoid.
        \param[in]	addG	Function pointer specifying the points that should be added to the goal set.
        \param[in]	addI	Function pointer specifying the points that should be added to the initial set.
    */
    template<class G_type, class I_type>
    void initializeReach(G_type addG, I_type addI) {

        TicToc timer;
        timer.tic();
        this->initializeSpec(&Gs_, addG);
        clog << "Gs_ initialized.\n";
        this->initializeSpec(&Is_, addI);
        clog << "Is_ initialized.\n";
        clog << "------------------------------------------------initializeSpec on Gs_, Is_, Os_: " << timer.toc() << " seconds.\n";

        // checking that specification is valid
        for (int i = 0; i < *this->system_->numAbs_; i++) {
            if ((Gs_[i]->symbolicSet_ & this->Os_[i]->symbolicSet_) != this->ddmgr_->bddZero()) {
                error("Error: G and O have nonzero intersection.\n");
            }
            if ((Is_[i]->symbolicSet_ & this->Os_[i]->symbolicSet_) != this->ddmgr_->bddZero()) {
                error("Error: I and O have nonzero intersection.\n");
            }
        }
        clog << "No obstacle problem with specification.\n";

        for (int i = 0; i < *this->system_->numAbs_; i++) {
            SymbolicSet* validZ = new SymbolicSet(*this->Xs_[i]);
            validZs_.push_back(validZ);
        }
        clog << "validZs_ initialized with empty domain.\n";

        for (int i = 0; i < *this->system_->numAbs_; i++) {
            SymbolicSet* validC = new SymbolicSet(*this->Xs_[i], *this->U_);
            validCs_.push_back(validC);
        }
        clog << "validCs_ initialized with empty domain.\n";

        saveVerify();
    }

    /*! Saves a snapshot of a controller and its domain into the sequence of final controllers and controller domains.
        \param[in] curAbs	0-index of the abstraction which the controller and controller domain that should be saved belong to.
    */
    void saveCZ(int curAbs) {
        SymbolicSet* C = new SymbolicSet(*this->Cs_[curAbs]);
        C->symbolicSet_ = this->Cs_[curAbs]->symbolicSet_;
        finalCs_.push_back(C);
        SymbolicSet* Z = new SymbolicSet(*this->Xs_[curAbs]);
        Z->symbolicSet_ = this->Zs_[curAbs]->symbolicSet_;
        finalZs_.push_back(Z);
    }

    /*! Saves and prints to log file some information related to the reachability/always-eventually specification. */
    void saveVerify() {
        printVec(Gs_, "G");
        printVec(Is_, "I");
        checkMakeDir("plotting");
        Gs_[*this->system_->numAbs_-1]->writeToFile("plotting/G.bdd");
        Is_[*this->system_->numAbs_-1]->writeToFile("plotting/I.bdd");
        checkMakeDir("G");
        saveVec(Gs_, "G/G");
    }
};

}

#endif /* UPFRONTREACH_HH_ */

