/*! \file Reach.hh
 *  Contains the Reach class. */

#ifndef REACH_HH_
#define REACH_HH_

#include "Adaptive.hh"

using namespace helper;

namespace scots {

/*! \class Reach
 *  \brief A class (derived from base Adaptive) that does adaptive multiscale abstraction-based synthesis for a reach specification.
 */
template<class X_type, class U_type>
class Reach: public virtual Adaptive<X_type, U_type> {
public:

    vector<SymbolicSet*> Gs_; /*!< Instance of *Xs_[i] containing goal states. */
    vector<SymbolicSet*> Is_; /*!< Instance of *Xs_[i] containing initial states. */
    vector<SymbolicSet*> validZs_; /*!< Contains winning states that act as savepoints. */
    vector<SymbolicSet*> validCs_; /*!< Controllers that act as savepoints. */
    vector<SymbolicSet*> finalCs_; /*!< Sequence of controllers that satisfy the specification. */
    vector<SymbolicSet*> finalZs_; /*!< Sequence of domains of finalCs_. */

    /*! Constructor for a Reach object. */
    Reach(int dimX, double* lbX, double* ubX, double* etaX, double tau,
          int dimU, double* lbU, double* ubU, double* etaU,
          double* etaRatio, double tauRatio, int nint,
          int numAbs, int readXX, int readAbs, char* logFile)
        : Adaptive<X_type, U_type>(dimX, lbX, ubX, etaX, tau,
                                   dimU, lbU, ubU, etaU,
                                   etaRatio, tauRatio, nint,
                                   numAbs, readXX, readAbs, logFile)

    {
    }

    /*! Deconstructor for a Reach object. */
    ~Reach() {
        deleteVec(Gs_);
        deleteVec(Is_);
        deleteVec(validZs_);
        deleteVec(validCs_);
        deleteVec(finalCs_);
        deleteVec(finalZs_);
    }

    /*! One iteration in an adaptive minimal fixed point.
        \param[in]		minToGoCoarser		Minimum number of iterations at an abstraction (not coarsest) before attemping to project to the next coarser abstraction.
        \param[in]		minToBeValid		Minimum number of iterations at an abstraction (not finest) after having gone coarser before a controller is declared valid and saved as a backup.
        \param[in]		earlyBreak			Whether the fixed point should end as soon as a strategy exists for the initial state.
        \param[in]		verbose				Whether messages detailing the steps of the algorithm should be printed to the log file.
        \param[in,out]	curAbs				0-index of the abstraction under consideration for the current iteration.
        \param[in,out]	iter				Counter for the total number of iterations of the minimal fixed point.
        \param[in,out]	justCoarsed			Status of the abstraction under consideration.
        \param[in,out]	iterCurAbs			Counter for the consecutive number of iterations for the current abstraction.
        \param[in,out]	reached				Whether the initial state has been declared winning.
        \param[in,out]	stop				Whether the fixed point should end.
    */
    void mu(int minToGoCoarser, int minToBeValid, int earlyBreak, int verbose, int* curAbs, int* iter, int* justCoarsed, int* iterCurAbs, int* reached, int* stop) {
        clog << "current abstraction: " << *curAbs << '\n';
        clog << "mu iteration: " << *iter << '\n';
        clog << "justCoarsed: " << *justCoarsed << '\n';
        clog << "iterCurAbs: " << *iterCurAbs << '\n';
        clog << "reached: " << *reached << '\n';
        clog << "controllers: " << finalCs_.size() << '\n';

        // get pre(Z)
        BDD preF = this->preC(this->Zs_[*curAbs]->symbolicSet_, this->Ts_[*curAbs]->symbolicSet_, *this->TTs_[*curAbs], *curAbs);
        // disjunct with goal (minimal fixed point)
        preF |= Gs_[*curAbs]->symbolicSet_;
        // find new {(x,u)}
        BDD N = preF & (!(this->Cs_[*curAbs]->symbolicSet_.ExistAbstract(*this->cubeU_)));
        // add new {(x,u)} to C
        this->Cs_[*curAbs]->symbolicSet_ |= N;
        // project onto Xs_[*curAbs]
        this->Zs_[*curAbs]->symbolicSet_ = preF.ExistAbstract(*this->notXvars_[*curAbs]);

        if (((this->Zs_[*curAbs]->symbolicSet_ & Is_[*curAbs]->symbolicSet_) != this->ddmgr_->bddZero()) && (*reached == 0)) {
            *reached = 1;
            if (earlyBreak == 1) {
                saveCZ(*curAbs);
                *stop = 1;
                clog << "\nTotal number of controllers: " << finalCs_.size() << '\n';
                return;
            }
        }

        if (*iter != 1) { // first iteration will always fail since we start with an empty Z
            if (N == this->ddmgr_->bddZero()) { // if iteration failed to yield new (x,u)
                if (*curAbs == this->numAbs_ - 1) { // if we're in the finest abstraction
                    if (*reached == 1) {
                        saveCZ(*curAbs);
                    }
                    else {
                    }
                    *stop = 1;
                    clog << "\nTotal number of controllers: " << finalCs_.size() << '\n';
                }
                else {
                    if (verbose) {
                        clog << "No new winning states; going finer.\n";
                    }
                    if (*justCoarsed == 1) { // current controller has not been declared valid
                        if (verbose) {
                            clog << "Current controller has not been declared valid.\n";
                            clog << "Removing last elements of finalCs_ and finalZs_.\n";
                            clog << "Resetting Zs_[" << *curAbs << "] and Cs_[" << *curAbs << "] from valids.\n";
                        }

                        // remove the last saved C and F from finalCs_ and Fs_
                        SymbolicSet* C = finalCs_.back();
                        finalCs_.pop_back();
                        delete(C);
                        SymbolicSet* Z = finalZs_.back();
                        finalZs_.pop_back();
                        delete(Z);
                        *justCoarsed = 0;

                        // reset the current winning states and controller
                        this->Zs_[*curAbs]->symbolicSet_ = validZs_[*curAbs]->symbolicSet_;
                        this->Cs_[*curAbs]->symbolicSet_ = validCs_[*curAbs]->symbolicSet_;

                    }
                    else {
                        if (verbose) {
                            clog << "Current controller has been declared valid.\n";
                            clog << "Saving Z and C of abstraction " << *curAbs << " into valids, finals.\n";
                            clog << "Finding inner approximation of Zs_[" << *curAbs << "] and copying into validZs_[" << *curAbs+1 << "].\n";
                        }
                        validZs_[*curAbs]->symbolicSet_ = this->Zs_[*curAbs]->symbolicSet_;
                        validCs_[*curAbs]->symbolicSet_ = this->Cs_[*curAbs]->symbolicSet_;

                        saveCZ(*curAbs);
                        this->innerFinerAligned(this->Zs_[*curAbs], this->Zs_[*curAbs+1], *curAbs);
                        validZs_[*curAbs+1]->symbolicSet_ = this->Zs_[*curAbs+1]->symbolicSet_;
                    }
                    *curAbs += 1;
                    *iterCurAbs = 1;
                }
            }
            else { // if there were new (x,u)
                if ((*justCoarsed == 1) && (*iterCurAbs >= minToBeValid)) {
                    if (verbose) {
                        clog << "Current controller now valid.\n";
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
                                clog << "Saving Z and C of abstraction " << *curAbs << " into validZs, validCs, finalZs, finalCs.\n";
                                clog << "Saving Z and C of abstraction " << *curAbs-1 << " into validZs_, validCs.\n";
                            }
                            saveCZ(*curAbs);
                            validZs_[*curAbs]->symbolicSet_ = this->Zs_[*curAbs]->symbolicSet_;
                            validCs_[*curAbs]->symbolicSet_ = this->Cs_[*curAbs]->symbolicSet_;
                            *justCoarsed = 1;
                            validZs_[*curAbs-1]->symbolicSet_ = this->Zs_[*curAbs-1]->symbolicSet_;
                            validCs_[*curAbs-1]->symbolicSet_ = this->Cs_[*curAbs-1]->symbolicSet_;
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
            if (verbose) {
                clog << "First iteration, nothing happens.\n";
            }
        }

        clog << "\n";
        *iter += 1;
    }

    /*!	Writes, should they exist, a sequence of controller and controller domain BDDs to directories 'C' and 'Z' respectively that satisfy the reachability specification.
        \param[in]  startAbs            0-index of the abstraction to start with.
        \param[in]	minToGoCoarser		Minimum number of growing fixed point iterations needed before an attempt to go to a coarser abstraction.
        \param[in]	minToBeValid		Minimum number of growing fixed point iterations needed before a controller is declared valid.
        \param[in]  earlyBreak          If 1, the synthesis ends as soon as I meets the domain of C.
        \param[in]	verbose				If 1, prints additional information during synthesis to the log file.

        \return     1 if controller(s) satisfying specification is/are synthesized; 0 otherwise.
    */
    int reach(int startAbs, int minToGoCoarser, int minToBeValid, int earlyBreak, int verbose = 1) {
        if ( (startAbs < 0) || (startAbs >= this->numAbs_) ) {
            error("Error: startAbs out of range.\n");
        }

        TicToc tt;
        tt.tic();

        // removing obstacles from transition relation
        this->removeFromTs(&(this->Os_));

        clog << "Os_ removed from Ts_, TTs_.\n";


        int curAbs = startAbs;

        int iter = 1;
        int justCoarsed = 0;
        int iterCurAbs = 1;
        int reached = 0;
        int stop = 0;

        while (1) {
            mu(minToGoCoarser, minToBeValid, earlyBreak, verbose, &curAbs, &iter, &justCoarsed, &iterCurAbs, &reached, &stop);
            if (stop) {
                break;
            }
        }
        if (reached) {
            clog << "Won.\n";

            checkMakeDir("C");
            saveVec(finalCs_, "C/C");
            checkMakeDir("Z");
            saveVec(finalZs_, "Z/Z");
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

    /*! Initializes objects specific to the following specifications: always-eventually, reach-while-avoid.
        \param[in]	addG	Function pointer specifying the points that should be added to the goal set.
        \param[in]	addI	Function pointer specifying the points that should be added to the initial set.
    */
    template<class G_type, class I_type>
    void initializeReach(G_type addG, I_type addI) {

        TicToc tt;
        tt.tic();
        this->initializeSpec(&Gs_, addG);
        clog << "Gs_ initialized.\n";
        this->initializeSpec(&Is_, addI);
        clog << "Is_ initialized.\n";
        clog << "------------------------------------------------initializeSpec on Gs_, Is_, Os_: ";
        tt.toc();

        // checking that specification is valid
        for (int i = 0; i < this->numAbs_; i++) {
            if ((Gs_[i]->symbolicSet_ & this->Os_[i]->symbolicSet_) != this->ddmgr_->bddZero()) {
                error("Error: G and O have nonzero intersection.\n");
            }
            if ((Is_[i]->symbolicSet_ & this->Os_[i]->symbolicSet_) != this->ddmgr_->bddZero()) {
                error("Error: I and O have nonzero intersection.\n");
            }
        }
        clog << "No obstacle problem with specification.\n";

        for (int i = 0; i < this->numAbs_; i++) {
            SymbolicSet* validZ = new SymbolicSet(*this->Xs_[i]);
            validZs_.push_back(validZ);
        }
        clog << "validZs_ initialized with empty domain.\n";

        for (int i = 0; i < this->numAbs_; i++) {
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
        Gs_[this->numAbs_-1]->writeToFile("plotting/G.bdd");
        Is_[this->numAbs_-1]->writeToFile("plotting/I.bdd");
        checkMakeDir("G");
        saveVec(Gs_, "G/G");

    }

    //    /*! Debugging function. Implementation of basic SCOTS reachability (single abstraction) in the Adaptive framework. */
    //    template<class G_type, class I_type, class O_type>
    //    void reachBasicDebug(int startAbs, G_type addG, I_type addI, O_type addO, int debug = 0) {

    //        int curAbs = startAbs;

    //        SymbolicSet G(*Xs_[curAbs]);
    //        addG(&G);
    //        SymbolicSet O(G);
    //        addO(&O);
    //        G.symbolicSet_ &= !(O.symbolicSet_);

    ////        G.printInfo(1);

    //        SymbolicSet F(*Xs_[curAbs], *U_);
    ////        F.symbolicSet_ = G.symbolicSet_ & U_->symbolicSet_;
    ////        F.printInfo(1);
    ////        U_->printInfo(1);

    //        SymbolicSet C(*Xs_[curAbs], *U_);

    //        BDD Z = ddmgr_->bddZero();

    //        SymbolicSet domain(C);
    //        domain.addGridPoints();

    //        int iter = 1;
    //        while ((C.symbolicSet_ & Is_[curAbs]->symbolicSet_) == ddmgr_->bddZero()) {

    //            //             project onto Xs_[curAbs]

    //            Z = F.symbolicSet_.ExistAbstract(*notXvars_[curAbs]);

    //            if (iter == debug) {
    //                SymbolicSet Zss(*Xs_[curAbs]);
    //                Zss.symbolicSet_ = Z;
    //                clog << "Z = F.symbolicSet_.ExistAbstract(notXvars)\n";
    //                Zss.printInfo(1);
    //            }

    //            // swap to X2s_[curAbs]
    //            Z = Z.Permute(permutesXtoX2_[curAbs]);

    //            SymbolicSet X2(*X2s_[curAbs]);
    //            X2.symbolicSet_ = Z;
    //            //            X2.printInfo(1);

    //            BDD nZ = !Z;

    //            if (iter == debug) {
    //                X2.symbolicSet_ = nZ;
    //                clog << "BDD nZ = !Z\n";
    //                X2.printInfo(1);
    //            }

    //            BDD T = Abs_[curAbs]->transitionRelation_;
    //            BDD Fbdd = T.AndAbstract(nZ, *notXUvars_[curAbs]);


    //            SymbolicSet Tss(C, X2);
    //            Tss.symbolicSet_ = T;
    //            //            Tss.printInfo(1);

    //            if (iter == debug) {
    //                SymbolicSet Fss(*Xs_[curAbs], *U_);
    //                Fss.symbolicSet_ = Fbdd;
    //                clog << "BDD Fbdd = T.AndAbstract(nZ, notXUvars)\n";
    //                Fss.printInfo(1);
    //            }

    //            BDD nF = !Fbdd;

    //            if (iter == debug) {
    //                SymbolicSet nFss(*Xs_[curAbs], *U_);
    //                nFss.symbolicSet_ = nF;
    //                clog << "BDD nF = !Fbdd\n";
    //                nFss.printInfo(1);
    //            }

    //            BDD TT = T.ExistAbstract(*cubesX2_[curAbs]);
    //            BDD preF = TT.AndAbstract(nF, *notXUvars_[curAbs]);

    //            if (iter == debug) {
    //                SymbolicSet preFss(*Xs_[curAbs], *U_);
    //                preFss.symbolicSet_ = preF;
    //                clog << "BDD preF = TT.AndAbstract(nF, notXUvars)\n";
    //                preFss.printInfo(1);
    //            }

    //            preF |= G.symbolicSet_;

    //            BDD N = preF & (!(C.symbolicSet_.ExistAbstract(*cubeU_)));

    //            if (iter == debug) {
    //                SymbolicSet Nss(*Xs_[curAbs], *U_);
    //                Nss.symbolicSet_ = N;
    //                clog << "BDD N = preF & (!(C.symbolicSet_.ExistAbstract(notXvars)))\n";
    //                Nss.printInfo(1);
    //            }

    //            C.symbolicSet_ |= N;

    //            F.symbolicSet_ = preF;

    //            iter++;

    //        }
    //        clog << "iterations: " << iter << '\n';
    //        C.printInfo(1);
    //        checkMakeDir("debug");
    //        C.writeToFile("debug/C.bdd");

    //        // normal SCOTS sanity check
    //        clog << "------------------------------normal SCOTS-----------------------------------\n";
    //        FixedPoint fp(Abs_[curAbs]);
    //        SymbolicSet C2(C);
    //        C2.symbolicSet_ = fp.reach(G.symbolicSet_, Is_[curAbs]->symbolicSet_, 1, 1);
    ////        C.symbolicSet_ &= domain.symbolicSet_;
    //        C2.printInfo(1);

    //        clog << "Cs are equal: " << ((C.symbolicSet_ <= C2.symbolicSet_) && (C2.symbolicSet_ <= C.symbolicSet_)) << "\n";
    ////        domain.printInfo(2);
    //    }

};

}

#endif /* REACH_HH_ */

