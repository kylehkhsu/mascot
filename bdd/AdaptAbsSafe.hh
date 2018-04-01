#ifndef ADAPTABSSAFE_HH_
#define ADAPTABSSAFE_HH_

#include <cstdio>
#include <vector>
#include <thread>

#include "SymbolicModelGrowthBound.hh"
#include "System.hh"
#include "Helper.hh"

using std::clog;
using std::freopen;
using std::string;
using std::vector;
using namespace helper;

namespace scots {

//enum SafeResult {CONVERGED, NOTCONVERGED};

class AdaptAbsSafe {
public:
    Cudd* ddmgr_; /*!< A single manager object common to all BDDs used in the program. */
    System* system_; /*!< Contains abstraction parameters. */
    vector<double*> etaXs_; /*!< *system_->numAbs_ x *system_->dimX_ matrix of state space grid spacings. */
    vector<double*> tau_; /*!< *system_->numAbs_ x 1 matrix of time steps. */
    vector<SymbolicSet*> Xs_; /*!< The *system_->numAbs_ "pre" state space abstractions, coarsest (0) to finest. */
	vector<SymbolicSet*> SafeInner_; /*!< Instance of *Xc_[i] containing inner approximation of safe states. */
	vector<SymbolicSet*> SafeOuter_; /*!< Instance of *Xc_[i] containing outer approximation of safe states. */
    vector<SymbolicSet*> Zs_; /*!< Instance of *Xs_[i] containing winning states. */
	vector<SymbolicSet*> Rs_; /*!< Instance of *Xs_[i] containing the states which were removed */
    vector<SymbolicSet*> X2s_; /*!< The *system_->numAbs_ "post" state space abstractions, coarsest (0) to finest. */
    SymbolicSet* U_; /*!< The single input space abstraction. */
    vector<SymbolicSet*> Cs_; /*!< Controller \subseteq *Xs_[i] x *U_. */
    int numBDDVars_; /*!< Total number of BDD variables used by ddmgr_. */
    vector<BDD*> cubesX_; /*!< *cubesX_[i] is a BDD with a single minterm of all 1s over the domain of *Xs_[i]. */
    vector<BDD*> cubesX2_; /*!< *cubesX2_[i] is a BDD with a single minterm of all 1s over the domain of *X2s_[i]. */
    BDD* cubeU_; /*!< A BDD with a single minterm of all 1s over the domain of the input SymbolicSet. */
    vector<BDD*> notXUvars_; /*!< notXUvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and *U_ and 1 otherwise. */
    vector<BDD*> notXvars_; /*!< notXvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and 1 otherwise. */
    vector<BDD*> cubesCoarser_;
    vector<SymbolicSet*> Ts_; /*!< Ts_[i] stores the transition relation of abstraction i. */
    vector<SymbolicSet*> TTs_; /*!< TTs_[i] is Ts_[i] with X2s_[i] existentially abtracted. */
    vector<int*> permutesXtoX2_; /*!< To transform a BDD over X variables to an equivalent one over X2 variables. */
    vector<int*> permutesX2toX_; /*!< To transform a BDD over X2 variables to an equivalent one over X variables. */
    vector<int*> permutesCoarser_;
    vector<int*> permutesFiner_;
    vector<OdeSolver*> solvers_; /*!< ODE solvers (Runge-Katta approximation) for each abstraction time step. */

    vector<SymbolicSet*> Gs_; /*!< Instance of *Xs_[i] containing goal states. */
    vector<SymbolicSet*> validZs_; /*!< Contains winning states that act as savepoints. */
    //vector<SymbolicSet*> validCs_; /*!< Controllers that act as savepoints. */
    vector<SymbolicSet*> finalCs_; /*!< Sequence of controllers that satisfy the specification. */
    vector<SymbolicSet*> finalZs_; /*!< Sequence of domains of finalCs_. */
    vector<int> finalAbs_; /*!< Sequence of abstraction layer indices of finalCs_. */
    vector<SymbolicSet*> Ds_; /*!< Instance of *Xs_[i] containing possible winning states. */
    vector<SymbolicSet*> innerDs_;
    vector<SymbolicSet*> computedDs_;
    vector<SymbolicSet*> savedZs_;
    
    /* Multi-threading for computation of the transition relation */
    int numPacket_; /*!< Number of packets/parallel threads for computation of the transition relation */
    vector<SymbolicSet*> splitDs_; /*!< Ds_[i] divided into packets, the length of the vector being same as the number of packets. */

    vector<SymbolicSet*> uTs_; /*!< Transition relations for exploring, i.e. all with coarsest space gridding but varying time sampling. */
	int recursion_; /*!< current recursion depth. */

    /*int minToGoCoarser_;
    int minToBeValid_;*/
    int earlyBreak_;
    int verbose_;

   /* int m_;
    int p_;*/

    double abTime_;
    double synTime_;

    /*!	Constructor for an AdaptAbsSafe object.
     *  \param[in]  logFile     Filename of program log.
     */
    AdaptAbsSafe(char* logFile) {
        freopen(logFile, "w", stderr);
        clog << logFile << '\n';

        abTime_ = 0;
        synTime_ = 0;
    }
    /*! Destructor for an AdaptAbsSafe object. */
    ~AdaptAbsSafe() {
        clog << "Abstraction construction time: " << abTime_ << " seconds.\n";
        clog << "Controller synthesis time: " << synTime_ << " seconds.\n";
        clog << "Abs + syn time: " << abTime_ + synTime_ << " seconds.\n";
        deleteVecArray(etaXs_);
        deleteVec(tau_);
        deleteVec(Xs_);
        deleteVec(SafeInner_);
		deleteVec(SafeOuter_);
        deleteVec(Zs_);
		deleteVec(Rs_);
        deleteVec(X2s_);
        delete U_;
        deleteVec(Cs_);
        deleteVec(cubesX_);
        deleteVec(cubesX2_);
        delete cubeU_;
        deleteVec(notXUvars_);
        deleteVec(notXvars_);
        deleteVec(cubesCoarser_);
        deleteVec(Ts_);
        deleteVec(TTs_);
        deleteVecArray(permutesXtoX2_);
        deleteVecArray(permutesX2toX_);
        deleteVecArray(permutesCoarser_);
        deleteVecArray(permutesFiner_);
        deleteVec(solvers_);
        deleteVec(Gs_);
        deleteVec(validZs_);
        //deleteVec(validCs_);
        deleteVec(finalCs_);
        deleteVec(finalZs_);
        deleteVec(Ds_);
        deleteVec(innerDs_);
        deleteVec(computedDs_);
        deleteVec(savedZs_);
        deleteVec(uTs_);
        fclose(stderr);
        delete ddmgr_;
    }

	void addNewLayer() {
		*system_->numAbs_ += 1;
	}

	template<class sys_type, class rad_type, class X_type, class U_type>
	void onTheFlySafe(sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        // Initialize SafeOuter_ with correct domain
        SafeOuter_[*system_->numAbs_ - 1]->symbolicSet_ = SafeInner_[*system_->numAbs_ - 1]->symbolicSet_; // outer and inner approximation of the safe set is same in the finest layer
        for (int i = *system_->numAbs_ - 1; i > 0; i--) { // compute the outer approximations of the safe sets accross all layers
            coarserOuter(SafeOuter_[i - 1], SafeOuter_[i], i - 1);
        }
        clog << "SafeOuter_ set according to specification.\n";
        
		// start by synthesizing all uTs_
//        for (int i = 0; i < *system_->numAbs_; i++){
//            Ds_[i]->addGridPoints();
//        }
        Ds_[0]->addGridPoints();
        
        // Debug purpose begin
        Ds_[0]->writeToFile("Ds0.bdd");
        Xs_[0]->writeToFile("Xs0.bdd");
        computedDs_[0]->writeToFile("computedDs_[0].bdd");
        SymbolicSet* D = new SymbolicSet(*computedDs_[0]);
        D->symbolicSet_ = !computedDs_[0]->symbolicSet_;
        D->writeToFile("nComputedDs_[0].bdd");
        // debug purpose end
		
		TicToc timer;
		timer.tic();
		computeExplorationAbstractions(sysNext, radNext, x, u);
		abTime_ += timer.toc();

		// Ts_[0] is uTs_[0]
		Ts_[0]->symbolicSet_ = uTs_[0]->symbolicSet_;
		TTs_[0]->symbolicSet_ = Ts_[0]->symbolicSet_.ExistAbstract(*notXUvars_[0]);

        // debug purpose
        cout << "\nTransition BDD of the coarsest layer:\n";
        Ts_[0]->printInfo(1);
        uTs_[0]->printInfo(1);
        //debug purpose end
        
		//// debug purpose
		//checkMakeDir("validZsInit");
		//saveVec(validZs_, "validZsInit/Zs_");

		// begin on-the-fly safety synthesis
		int print = 1; // print progress on command line
        int recursion = 0;
		onTheFlySafeRecurse(sysNext, radNext, x, u, recursion, print);

		// to ensure that the controllers for the coarser layers admit control inputs that allow visit to finer layer states. Note that this doesn't change the controller domain.
		for (int ab = 0; ab < *system_->numAbs_; ab++) {
			BDD cPreZ = cPre(validZs_[ab]->symbolicSet_, ab);
			Cs_[ab]->symbolicSet_ = cPreZ & SafeInner_[ab]->symbolicSet_;
		}
		saveCZ();

		int NoController = 1;
		for (int ab = 0; ab < *system_->numAbs_; ab++) {
			if (Cs_[ab]->symbolicSet_ != ddmgr_->bddZero()) {
				NoController = 0;
				// clog << "controllers: " << finalCs_.size() << '\n';
				checkMakeDir("C");
				saveVec(finalCs_, "C/C");
				checkMakeDir("Z");
				saveVec(finalZs_, "Z/Z");
			}
		}

		if (NoController == 1) {
			clog << "Empty controller domain.\n";
		}
		
        // debug purpose
        cout << "\nTransition BDD of the coarsest layer:\n";
        Ts_[0]->printInfo(1);
//        cout << "\nTransition BDD of the second coarsest layer:\n";
//        Ts_[1]->printInfo(1);
        // debug purpose end
		checkMakeDir("T");
		saveVec(Ts_, "T/T");
		clog << "Wrote Ts_ to file.\n";
//        return;
	}

	template<class sys_type, class rad_type, class X_type, class U_type>
	void onTheFlySafeNoAux(sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        TicToc timer;
        timer.tic();

		// Initialize SafeOuter_ with correct domain
		SafeOuter_[*system_->numAbs_ - 1]->symbolicSet_ = SafeInner_[*system_->numAbs_ - 1]->symbolicSet_; // outer and inner approximation of the safe set is same in the finest layer
		for (int i = *system_->numAbs_ - 1; i > 0; i--) { // compute the outer approximations of the safe sets accross all layers
			coarserOuter(SafeOuter_[i - 1], SafeOuter_[i], i - 1);
		}
		clog << "SafeOuter_ set according to specification.\n";

		// start by synthesizing all uTs_
		//        for (int i = 0; i < *system_->numAbs_; i++){
		//            Ds_[i]->addGridPoints();
		//        }
		Ds_[0]->addGridPoints();		
	
		SymbolicModelGrowthBound<X_type, U_type> abstraction(Ds_[0], U_, X2s_[0]); // coarsest state gridding
		abstraction.computeTransitionRelation(sysNext, radNext, *solvers_[0]); // use solver with time step corresponding to that of each layer
		uTs_[0]->symbolicSet_ = abstraction.transitionRelation_; // add to transition relation
		
		// Ts_[0] is uTs_[0]
		Ts_[0]->symbolicSet_ = uTs_[0]->symbolicSet_;
		TTs_[0]->symbolicSet_ = Ts_[0]->symbolicSet_.ExistAbstract(*notXUvars_[0]);

        abTime_ += timer.toc();
		//// debug purpose
		//checkMakeDir("validZsInit");
		//saveVec(validZs_, "validZsInit/Zs_");

		// begin on-the-fly safety synthesis
		int print = 1; // print progress on command line
		int recursion = 0;
        checkMakeDir("K");
		onTheFlySafeRecurseNoAux(sysNext, radNext, x, u, recursion, print);

		// to ensure that the controllers for the coarser layers admit control inputs that allow visit to finer layer states. Note that this doesn't change the controller domain.

        timer.tic();
		for (int ab = 0; ab < *system_->numAbs_; ab++) {
			BDD cPreZ = cPre(validZs_[ab]->symbolicSet_, ab);
			Cs_[ab]->symbolicSet_ = cPreZ & SafeInner_[ab]->symbolicSet_;
		}
        synTime_ += timer.toc();
		saveCZ();

		int NoController = 1;
		for (int ab = 0; ab < *system_->numAbs_; ab++) {
			if (Cs_[ab]->symbolicSet_ != ddmgr_->bddZero()) {
				NoController = 0;
				// clog << "controllers: " << finalCs_.size() << '\n';
				checkMakeDir("C");
				saveVec(finalCs_, "C/C");
				checkMakeDir("Z");
				saveVec(finalZs_, "Z/Z");
			}
		}

		if (NoController == 1) {
			clog << "Empty controller domain.\n";
		}

		checkMakeDir("T");
		saveVec(Ts_, "T/T");
		clog << "Wrote Ts_ to file.\n";
		//        return;
	}
    
	template<class sys_type, class rad_type, class X_type, class U_type>
	void computeExplorationAbstractions(sys_type sysNext, rad_type radNext, X_type x, U_type u) {
		for (int ab = 0; ab < *system_->numAbs_; ab++) {
			SymbolicModelGrowthBound<X_type, U_type> abstraction(Ds_[0], U_, X2s_[0]); // coarsest state gridding
			abstraction.computeTransitionRelation(sysNext, radNext, *solvers_[ab]); // use solver with time step corresponding to that of each layer
			uTs_[ab]->symbolicSet_ = abstraction.transitionRelation_; // add to transition relation
		}
	}

	template<class sys_type, class rad_type, class X_type, class U_type>
	void onTheFlySafeRecurse(sys_type sysNext, rad_type radNext, X_type x, U_type u, int recursion, int print = 0) {
		recursion += 1;
		bool CONVERGED;
		clog << '\n';
		clog << "current recursion depth: " << recursion << '\n';

		if (print) {
			cout << '\n';
			cout << "current recursion depth: " << recursion << '\n';
		}

		for (int ab = 0; ab < *system_->numAbs_; ab++) {
            clog << "\t current abstraction: " << ab << '\n';
            if (print)
                cout << "\t current abstraction: " << ab << '\n';
			if (ab > 0) {
                TicToc timer;
                timer.tic();
				ExpandAbstraction(ab, sysNext, radNext, x, u);
                abTime_ += timer.toc();
			}
            TicToc timer;
            timer.tic();
			safeOne(ab, print);
            synTime_ += timer.toc();
		}


        for (int ab = 1; ab <= *system_->numAbs_ - 1; ab++) {
            BDD X = Zs_[ab]->symbolicSet_;
            finer(Zs_[ab - 1], Zs_[ab], ab - 1);
            Zs_[ab]->symbolicSet_ |= X;
        }

		// Debug purpose
		/*if (recursion == 1) {
			checkMakeDir("Zs1");
			saveVec(Zs_, "Zs1/Zs_");
			checkMakeDir("validZs1");
			saveVec(validZs_, "validZs1/Zs_");
		}
		else if (recursion == 2) {
			checkMakeDir("Zs2");
			saveVec(Zs_, "Zs2/Zs_");
			checkMakeDir("validZs2");
			saveVec(validZs_, "validZs2/Zs_");
		}
		else if (recursion == 3) {
			checkMakeDir("Zs3");
			saveVec(Zs_, "Zs3/Zs_");
			checkMakeDir("validZs3");
			saveVec(validZs_, "validZs3/Zs_");
		}
		else if (recursion == 4) {
			checkMakeDir("Zs4");
			saveVec(Zs_, "Zs4/Zs_");
			checkMakeDir("validZs4");
			saveVec(validZs_, "validZs4/Zs_");
		}*/
		

		//Convergence criteria: validZs_ of the finest layer is same as Zs_ of the finest layer
		if (validZs_[*system_->numAbs_ - 1]->symbolicSet_ == Zs_[*system_->numAbs_ - 1]->symbolicSet_) {
			CONVERGED = 1;
		} 
		else {
			CONVERGED = 0;
		}

		

		for (int ab = *system_->numAbs_ - 1; ab >= 0; ab--) {
			validZs_[ab]->symbolicSet_ = ddmgr_->bddZero(); // just to make sure there is no interference with previous result
		}

		validZs_[*system_->numAbs_ - 1]->symbolicSet_ = Zs_[*system_->numAbs_ - 1]->symbolicSet_; // save the winning states for Cpre in next iteration
		for (int ab = *system_->numAbs_ - 1; ab > 0; ab--) {
			coarserInner(validZs_[ab-1], validZs_[ab], ab-1);
		}

		if (CONVERGED == 1) {
			return;
		} 
		else {
            for (int ab = 0; ab <= *system_->numAbs_ - 1; ab++){
                Zs_[ab]->symbolicSet_ = ddmgr_->bddZero();
            }
			onTheFlySafeRecurse(sysNext, radNext, x, u, recursion, print);
			return;
		}
	}

	template<class sys_type, class rad_type, class X_type, class U_type>
	void onTheFlySafeRecurseNoAux(sys_type sysNext, rad_type radNext, X_type x, U_type u, int recursion, int print = 0) {
        TicToc timer;
		recursion += 1;
		bool CONVERGED;

        if (print) {
		    clog << '\n';
		    clog << "current recursion depth: " << recursion << '\n';

		
			cout << '\n';
			cout << "current recursion depth: " << recursion << '\n';
		}

		/*Ds_[0]->symbolicSet_ = validZs_[0]->symbolicSet_;*/
		for (int ab = 0; ab < *system_->numAbs_; ab++) {			
			if (print) {
                clog << "\t current abstraction: " << ab << '\n';
				cout << "\t current abstraction: " << ab << '\n';
            }
			if (ab > 0) {
				//Ds_[ab-1]->symbolicSet_ = !Zs_[ab - 1]->symbolicSet_ & validZs_[ab - 1]->symbolicSet_;
			    timer.tic();
				Ds_[*system_->numAbs_ - 1]->symbolicSet_ = validZs_[*system_->numAbs_ - 1]->symbolicSet_;
				for (int ii = *system_->numAbs_ - 1; ii >= ab; ii--) {
					coarserOuter(Ds_[ii - 1], Ds_[ii], ii - 1);
				}
				Ds_[ab-1]->symbolicSet_ &= Rs_[ab-1]->symbolicSet_;
				finer(Ds_[ab - 1], Ds_[ab], ab - 1);
				computeAbstraction(ab, sysNext, radNext, x, u);
				abTime_ += timer.toc();
			}
			timer.tic();
			safeOne(ab, print);
			synTime_ += timer.toc();

			/*if (ab > 0) {
				BDD X = Zs_[ab]->symbolicSet_;
				finer(Zs_[ab - 1], Zs_[ab], ab - 1);
				Zs_[ab]->symbolicSet_ |= X;
			}*/
		}

		timer.tic();
		for (int ab = 1; ab <= *system_->numAbs_ - 1; ab++) {
			BDD X = Zs_[ab]->symbolicSet_;
			finer(Zs_[ab - 1], Zs_[ab], ab - 1);
			Zs_[ab]->symbolicSet_ |= X;
		}
        synTime_ += timer.toc();

		// Debug purpose
		/*if (recursion == 1) {
		checkMakeDir("Zs1");
		saveVec(Zs_, "Zs1/Zs_");
		checkMakeDir("validZs1");
		saveVec(validZs_, "validZs1/Zs_");
		}
		else if (recursion == 2) {
		checkMakeDir("Zs2");
		saveVec(Zs_, "Zs2/Zs_");
		checkMakeDir("validZs2");
		saveVec(validZs_, "validZs2/Zs_");
		}
		else if (recursion == 3) {
		checkMakeDir("Zs3");
		saveVec(Zs_, "Zs3/Zs_");
		checkMakeDir("validZs3");
		saveVec(validZs_, "validZs3/Zs_");
		}
		else if (recursion == 4) {
		checkMakeDir("Zs4");
		saveVec(Zs_, "Zs4/Zs_");
		checkMakeDir("validZs4");
		saveVec(validZs_, "validZs4/Zs_");
		}*/

		timer.tic();
		//Convergence criteria: validZs_ of the finest layer is same as Zs_ of the finest layer
		if (validZs_[*system_->numAbs_ - 1]->symbolicSet_ == Zs_[*system_->numAbs_ - 1]->symbolicSet_) {
			CONVERGED = 1;
		}
		else {
			CONVERGED = 0;
		}
		for (int ab = *system_->numAbs_ - 1; ab >= 0; ab--) {
			validZs_[ab]->symbolicSet_ = ddmgr_->bddZero(); // just to make sure there is no interference with previous result
		}

		validZs_[*system_->numAbs_ - 1]->symbolicSet_ = Zs_[*system_->numAbs_ - 1]->symbolicSet_; // save the winning states for Cpre in next iteration
		for (int ab = *system_->numAbs_ - 1; ab > 0; ab--) {
			coarserInner(validZs_[ab - 1], validZs_[ab], ab - 1);
		}
        synTime_ += timer.toc();
        
        // debug purpose
        string Str = "K/K";
        Str += std::to_string(recursion);
        checkMakeDir(Str);
        saveVec(validZs_, Str);
        // debug purpose end

		if (CONVERGED == 1) {
			return;
		}
		else {
		    timer.tic();
			for (int ab = 0; ab <= *system_->numAbs_ - 1; ab++) {
				Zs_[ab]->symbolicSet_ = ddmgr_->bddZero();
				/*Ds_[ab]->symbolicSet_ = validZs_[ab]->symbolicSet_;*/
			}
            synTime_ += timer.toc();
			onTheFlySafeRecurseNoAux(sysNext, radNext, x, u, recursion, print);
			return;
		}
	}

	template<class sys_type, class rad_type, class X_type, class U_type>
	void ExpandAbstraction(int nextAb, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
		Ds_[*system_->numAbs_ - 1]->symbolicSet_ = validZs_[*system_->numAbs_ - 1]->symbolicSet_; // save the winning states for Pre in the next iteration
		for (int ab = *system_->numAbs_ - 1; ab > 0; ab--) {
			coarserOuter(Ds_[ab - 1], Ds_[ab], ab - 1);
		}
		BDD uPreD = uPre(Ds_[0]->symbolicSet_, nextAb);
		Ds_[0]->symbolicSet_ = uPreD & SafeOuter_[0]->symbolicSet_;
		for (int ab = 0; ab < nextAb - 1; ab++) {
			finer(Ds_[ab], Ds_[ab + 1], ab);
		}
		Ds_[nextAb - 1]->symbolicSet_ &= Rs_[nextAb - 1]->symbolicSet_;
		finer(Ds_[nextAb - 1], Ds_[nextAb], nextAb - 1);
        computeAbstraction(nextAb, sysNext, radNext, x, u);
//        computeAbstractionInParallel(nextAb, sysNext, radNext, x, u);
	}
    
    // one step safety fixed-point
    void safeOne(int ab, int print = 0) {
        int i = 1;
        
        BDD cPreZ = cPre(validZs_[ab]->symbolicSet_, ab);
        BDD C = cPreZ & SafeInner_[ab]->symbolicSet_;
        if (ab == 0) {
            Rs_[ab]->symbolicSet_ = (!(C.ExistAbstract(*cubeU_))) & validZs_[ab]->symbolicSet_; // the states which were removed (in the present iteration only)
        }
        else {
            finer(Rs_[ab - 1], Rs_[ab], ab - 1);
            Rs_[ab]->symbolicSet_ &= (!(C.ExistAbstract(*cubeU_))); // & Rs_[ab - 1]->symbolicSet_; // the states which were removed (in the present iteration only)
        }
        Cs_[ab]->symbolicSet_ = C;
        Zs_[ab]->symbolicSet_ |= Cs_[ab]->symbolicSet_.ExistAbstract(*cubeU_);
    }

//    // one step safety fixed-point
//    void safeOne(int ab, int print = 0) {
//        int i = 1;
//
//        BDD cPreZ = cPre(validZs_[ab]->symbolicSet_, ab);
//        BDD C = cPreZ.ExistAbstract(*cubeU_) & SafeInner_[ab]->symbolicSet_;
//        if (ab == 0) {
//            Rs_[ab]->symbolicSet_ = !C & validZs_[ab]->symbolicSet_; // the states which were removed (in the present iteration only)
//        }
//        else {
//            finer(Rs_[ab - 1], Rs_[ab], ab - 1);
//            Rs_[ab]->symbolicSet_ &= !C; // & Rs_[ab - 1]->symbolicSet_; // the states which were removed (in the present iteration only)
//        }
//        Cs_[ab]->symbolicSet_ &= !Rs_[ab]->symbolicSet_;
//        Zs_[ab]->symbolicSet_ |= Cs_[ab]->symbolicSet_.ExistAbstract(*cubeU_);
//    }

    // testing computeAbstraction
    void test() {
        Ds_[0]->addGridPoints();
    }

    // abstraction synthesis will only iterate over the elements in the domain of the state space SymbolicSet (first argument in constructing abstraction)
    template<class sys_type, class rad_type, class X_type, class U_type>
    void computeAbstraction(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u, int print = 0) {
        BDD D = Ds_[ab]->symbolicSet_ & (!computedDs_[ab]->symbolicSet_); // don't repeat sampling
        computedDs_[ab]->symbolicSet_ |= Ds_[ab]->symbolicSet_; // update computed part of transition relation
        Ds_[ab]->symbolicSet_ = D;

        if (ab == 1) {
            checkMakeDir("D1");
            if (!(std::ifstream("D1/1.bdd"))) {
                Ds_[1]->writeToFile("D1/1.bdd");
            } else if (!(std::ifstream("D1/2.bdd"))) {
                Ds_[1]->writeToFile("D1/2.bdd");
            } else if (!(std::ifstream("D1/3.bdd"))) {
                Ds_[1]->writeToFile("D1/3.bdd");
            } else if (!(std::ifstream("D1/4.bdd"))) {
                Ds_[1]->writeToFile("D1/4.bdd");
            } else if (!(std::ifstream("D1/5.bdd"))) {
                Ds_[1]->writeToFile("D1/5.bdd");
            } else if (!(std::ifstream("D1/6.bdd"))) {
                Ds_[1]->writeToFile("D1/6.bdd");
            }
        }
        if (ab == 2) {
            checkMakeDir("D2");
            if (!(std::ifstream("D2/1.bdd"))) {
                Ds_[2]->writeToFile("D2/1.bdd");
            } else if (!(std::ifstream("D2/2.bdd"))) {
                Ds_[2]->writeToFile("D2/2.bdd");
            } else if (!(std::ifstream("D2/3.bdd"))) {
                Ds_[2]->writeToFile("D2/3.bdd");
            } else if (!(std::ifstream("D2/4.bdd"))) {
                Ds_[2]->writeToFile("D2/4.bdd");
            } else if (!(std::ifstream("D2/5.bdd"))) {
                Ds_[2]->writeToFile("D2/5.bdd");
            } else if (!(std::ifstream("D2/6.bdd"))) {
                Ds_[2]->writeToFile("D2/6.bdd");
            }
        }

        SymbolicModelGrowthBound<X_type, U_type> abstraction(Ds_[ab], U_, X2s_[ab]);
        abstraction.computeTransitionRelation(sysNext, radNext, *solvers_[ab], 0); // was hard-coded to *solvers_[0], source of "tunneling" bug
        if (print)
		    (abstraction.getTransitionRelation()).printInfo(1);
        if (abstraction.transitionRelation_ != ddmgr_->bddZero()) { // no point adding/displaying if nothing was added
            Ts_[ab]->symbolicSet_ |= abstraction.transitionRelation_; // add to transition relation
//            BDD O2 = Os_[ab]->symbolicSet_.Permute(permutesXtoX2_[ab]); // causes unsound behavior, leaving here as warning
//            Ts_[ab]->symbolicSet_ &= !O2;
//            Ts_[ab]->printInfo(1);
            TTs_[ab]->symbolicSet_ = Ts_[ab]->symbolicSet_.ExistAbstract(*notXUvars_[ab]);
        }
    }

    //template<class sys_type, class rad_type, class X_type, class U_type>
    //void computeAbstractionInParallel(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
    //    BDD D = Ds_[ab]->symbolicSet_ & (!computedDs_[ab]->symbolicSet_); // don't repeat sampling
    //    computedDs_[ab]->symbolicSet_ |= Ds_[ab]->symbolicSet_; // update computed part of transition relation
    //    Ds_[ab]->symbolicSet_ = D;
    //    
    //    initializeParallel(ab);
    //    vector<SymbolicModelGrowthBound<X_type,U_type>> splitAbstraction;
    //    Ds_[ab]->splitGridPoints(splitDs_, numPacket_);
    //    for (int packet = 0; packet < numPacket_; packet++) {
    //        SymbolicModelGrowthBound<X_type, U_type> abs(splitDs_[packet], U_, X2s_[ab]);
    //        splitAbstraction.push_back(abs);
    //    }
    //    std::thread t[numPacket_];
    //    for (int packet = 0; packet < numPacket_; packet++) {
    //        t[packet] = std::thread(this->abstractionFromThread(splitAbstraction, sysNext, radNext, ab, packet));
    //    }
    //    
    //    for (int packet = 0; packet < numPacket_; packet++) {
    //        t[packet].join();
    //    }
    //    
    //    SymbolicModelGrowthBound<X_type, U_type> abstraction(Ds_[ab], U_, X2s_[ab]);
    //    for (int packet = 0; packet < numPacket_; packet++) {
    //        abstraction.transitionRelation_ |= splitAbstraction[packet].transitionRelation_;
    //    }
    //    deleteVec(splitDs_);
    //    
    //    if (abstraction.transitionRelation_ != ddmgr_->bddZero()) { // no point adding/displaying if nothing was added
    //        Ts_[ab]->symbolicSet_ |= abstraction.transitionRelation_; // add to transition relation
    //        //            BDD O2 = Os_[ab]->symbolicSet_.Permute(permutesXtoX2_[ab]); // causes unsound behavior, leaving here as warning
    //        //            Ts_[ab]->symbolicSet_ &= !O2;
    //        //            Ts_[ab]->printInfo(1);
    //        TTs_[ab]->symbolicSet_ = Ts_[ab]->symbolicSet_.ExistAbstract(*notXUvars_[ab]);
    //    }
    //}
    //
    //template<class SymbolicModelGrowthBound, class sys_type, class rad_type>
    //void abstractionFromThread(vector<SymbolicModelGrowthBound*> splitAbstraction, sys_type sysNext, rad_type radNext, int ab, int packet) {
    //    splitAbstraction[packet].computeTransitionRelation(sysNext, radNext, *solvers_[ab], 1);
    //}
    
//    /*! Divide Ds_[ab] into same-sized packets */
//    void SplitInPakets(int ab){
//        initializeParallel(ab);
//        for (Ds_[ab]->begin(); !Ds_[ab]->done(); Ds_[ab]->next()){
//            for (int packet = 0; packet < numPacket_; packet++) {
//                splitDs_[packet]->symbolicSet_ |= Ds_[ab]->currentMinterm();
//            }
//        }
////        deleteVec(splitDs_);
//    }
    

    /*! Calculates the cooperative predecessor of a given set with respect to the un-restricted transition relation at the coarsest level of abstraction.
     *  \param[in]  Z       The given set.
     *  \return     BDD containing {x} for which there exists a u s.t. there exists a post state of (x,u) in Z.
     */
    BDD uPre(BDD Z, int ab) {
        BDD Z2 = Z.Permute(permutesXtoX2_[0]);
        BDD uPreZ = uTs_[ab]->symbolicSet_.AndAbstract(Z2, *notXUvars_[0]);
        uPreZ = uPreZ.ExistAbstract(*notXvars_[0]);
        return uPreZ;
    }


    /*! Calculates the enforceable predecessor of the given set with respect to the transition relation at the specified level of abstraction.
     *  \param[in]  Z           The winning set.
     *  \param[in]  curAbs      0-index of the current abstraction.
     *  \return     BDD containing {(x,u)} for which all post states are in Z.
     */
    BDD cPre(BDD Z, int curAbs) {
        // swap to X2s_[curAbs]
        BDD Z2 = Z.Permute(permutesXtoX2_[curAbs]);
        // posts outside of winning set
        BDD nZ2 = !Z2;
        // {(x,u)} with some posts outside of winning set
        BDD Fbdd = Ts_[curAbs]->symbolicSet_.AndAbstract(nZ2, *notXUvars_[curAbs]);
        // {(x,u)} with no posts outside of winning set
        BDD nF = !Fbdd;
        // get rid of junk
        BDD preF = TTs_[curAbs]->symbolicSet_.AndAbstract(nF, *notXUvars_[curAbs]);        
        return preF;
    }

    void finer(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
//        std::cout << "\n nextAb: " << c;
//        std::cout << "\n State space of Ds_[0]: \n" << Zc->symbolicSet_;
//        std::cout << "\n State space of Ds_[1]: \n" << Zf->symbolicSet_;
        Zf->addGridPoints();
//        std::cout << "\n State space of Ds_[1] (after addGP): \n" << Zf->symbolicSet_;
        Zf->symbolicSet_ &= Zc->symbolicSet_.Permute(permutesFiner_[c]);
//        std::cout << "\n State space of Ds_[1] (after conj): \n" << Zf->symbolicSet_;
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

    /*! Initializes data members.
     *  \param[in]  system      Contains abstraction parameters.
     *  \param[in]	addS	Function pointer specifying the points that should be added to the safe set.
     *  \param[in]  num     Number of parallel threads
     */
    template<class S_type>
    void initialize(System* system, S_type addS, int num=1) {
        ddmgr_ = new Cudd;
        system_ = system;
        numPacket_ = num;

        TicToc timer;
        timer.tic();

        initializeEtaTau();
        initializeSolvers();
        initializeSymbolicSets(addS);
        initializeNumBDDVars();
        initializePermutes();
        initializeCubes();
        initializeNotVars();
        initializeProjs();
        verifySave();

        abTime_ += timer.toc();
    }

    /*! Initializes SymbolicSet data members.
     *  \param[in]	addS	Function pointer specifying the points that should be added to the inner approx of the safe set.
     */
    template<class S_type>
    void initializeSymbolicSets(S_type addS) {
        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* X = new SymbolicSet(*ddmgr_, *system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[i], tau_[i][0]);
            X->addGridPoints();
            Xs_.push_back(X);
        }
        clog << "Xs_ initialized with full domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* X2 = new SymbolicSet(*Xs_[i], 1);
            X2->addGridPoints();
            X2s_.push_back(X2);
        }
        clog << "X2s_ initialized with full domain.\n";

        U_ = new SymbolicSet(*ddmgr_, *system_->dimU_, system_->lbU_, system_->ubU_, system_->etaU_, 0);
        U_->addGridPoints();
        clog << "U_ initialized with full domain.\n";

        initializeSpec(&SafeInner_, addS);
        clog << "SafeInner_ initialized according to specification.\n";
		
		for (int i = 0; i < *system_->numAbs_; i++) {
			SymbolicSet* SafeOuter = new SymbolicSet(*Xs_[i]);
			SafeOuter_.push_back(SafeOuter);
		}
        clog << "SafeOuter_ initialized with empty domain, and will be set with correct domain later.";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* Z = new SymbolicSet(*Xs_[i]);
			//Z->addGridPoints();
            Zs_.push_back(Z);

        }
        clog << "Zs_ initialized with empty domain.\n";

		for (int i = 0; i < *system_->numAbs_; i++) {
			SymbolicSet* R = new SymbolicSet(*Xs_[i]);
			Rs_.push_back(R);

		}
		clog << "Rs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* C = new SymbolicSet(*Xs_[i], *U_);
			C->addGridPoints();
            Cs_.push_back(C);
        }
        clog << "Cs_ initialized with full domain.\n";

        
        // checking that specification is valid
        if (SafeInner_[*system_->numAbs_ - 1]->symbolicSet_ == ddmgr_->bddZero()) {
			 error("Error: Safe set is empty.\n");
            clog << "Safe set is empty.\n";
        }

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* validZ = new SymbolicSet(*Xs_[i]);
			validZ->addGridPoints();
            validZs_.push_back(validZ);
        }
        clog << "validZs_ initialized with full domain.\n";

        /*for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* validC = new SymbolicSet(*Xs_[i], *U_);
			validC->addGridPoints();
            validCs_.push_back(validC);
        }
        clog << "validCs_ initialized with full domain.\n";*/

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* D = new SymbolicSet(*Xs_[i]);
            Ds_.push_back(D);
        }
        clog << "Ds_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* innerD = new SymbolicSet(*Xs_[i]);
            innerDs_.push_back(innerD);
        }
        clog << "innerDs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* computedD = new SymbolicSet(*Xs_[i]);
//            computedD->addGridPoints();
//            computedD->symbolicSet_ = ddmgr_->bddZero();
            computedDs_.push_back(computedD);
        }
        clog << "computedDs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* savedZ = new SymbolicSet(*Xs_[i]);
            savedZs_.push_back(savedZ);
        }
        clog << "savedZs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* T = new SymbolicSet(*Cs_[i], *X2s_[i]);
//            T->symbolicSet_ = ddmgr_->bddZero();
            Ts_.push_back(T);
        }
        clog << "Ts_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* TT = new SymbolicSet(*Cs_[i]);
//            TT->symbolicSet_ = ddmgr_->bddZero();
            TTs_.push_back(TT);
        }
        clog << "TTs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* uT = new SymbolicSet(*Ts_[0]); // exploration transition relations are all with coarsest gridding
            uT->symbolicSet_ = ddmgr_->bddZero();
            uTs_.push_back(uT);
        }
        clog << "uT_ initialized with empty domain.\n";
    }

    /*! Initializes a vector of SymbolicSets that is an instance of Xs_ containing points as specified by addSpec.
     *  \param[in] vec      Vector of SymbolicSets to initialize.
     *  \param[in] addSpec  Function pointer specifying the points to add to each element of vec.
     */
    template<class vec_type, class spec_type>
    void initializeSpec(vec_type* vec, spec_type addSpec) {
        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* Xinstance = new SymbolicSet(*Xs_[i]);
            addSpec(Xinstance);
            vec->push_back(Xinstance);
        }
    }

    /*! Initializes the BDDs useful for existential abstraction. Must be called only after initializing Xs, U, and X2s. */
    void initializeNotVars() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            BDD* notXUvars = new BDD;
            BDD* notXvars = new BDD;
            *notXUvars = ddmgr_->bddOne();
            for (int j = 0; j < *system_->numAbs_; j++) {
                *notXUvars &= *cubesX2_[j];
                if (i != j) {
                    *notXUvars &= *cubesX_[j];
                }
            }
            *notXvars = *notXUvars & *cubeU_;

            notXUvars_.push_back(notXUvars);
            notXvars_.push_back(notXvars);
        }

        clog << "Initialized notVars.\n";
    }

    /*! Initializes the abstractions' state space grid parameters and time sampling parameters. */
    void initializeEtaTau() {
        for (int i = 0; i < *system_->dimX_; i++) {
            if (system_->etaRatio_[i] != 1 && system_->etaRatio_[i] != 2) {
                error("Error: unsupported etaRatio, which must consist of only 1s or 2s.\n");
            }
        }

        double* etaCur = new double[*system_->dimX_];
        for (int i = 0; i < *system_->dimX_; i++) {
            etaCur[i] = system_->etaX_[i];
        }
        double* tauCur = new double;
        *tauCur = *system_->tau_;

        for (int i = 0; i < *system_->numAbs_; i++) {
            double* etaX = new double[*system_->dimX_];
            double* taui = new double;
            for (int j = 0; j < *system_->dimX_; j++) {
                etaX[j] = etaCur[j];
            }
            *taui = *tauCur;
            etaXs_.push_back(etaX);
            tau_.push_back(taui);

            for (int j = 0; j < *system_->dimX_; j++) {
                etaCur[j] /= system_->etaRatio_[j];
            }
            *tauCur /= *system_->tauRatio_;
        }

        clog << "Initialized etaXs_, tau_.\n";

        delete[] etaCur;
        delete tauCur;
    }

    /*! Initializes the Runge-Katta ODE solvers.
     */
    void initializeSolvers() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            OdeSolver* solver = new OdeSolver(*system_->dimX_, *system_->nSubInt_, *tau_[i]);
            solvers_.push_back(solver);
        }
        clog << "Initialized solvers.\n";
    }

    /*! Initializes the number of total distinct BDD variables in use. */
    void initializeNumBDDVars() {
        numBDDVars_ = 0;
        for (int i = 0; i < *system_->numAbs_; i++) {
            numBDDVars_ += Xs_[i]->nvars_;
            numBDDVars_ += X2s_[i]->nvars_;
        }
        numBDDVars_ += U_->nvars_;

        clog << "Number of BDD variables: " << numBDDVars_ << '\n';
    }

    void initializeProjs() {
        int ones = 0;
        for (int i = 0; i < *system_->dimX_; i++) {
            if (system_->etaRatio_[i] == 1) {
                ones++;
            }
        }

        for (int i = 1; i < *system_->numAbs_; i++) {
            BDD* varsCoarser = new BDD[*system_->dimX_ - ones];
            int ind = 0;
            for (int j = 0; j < *system_->dimX_; j++) {
                if (system_->etaRatio_[j] == 2) {
                    varsCoarser[ind] = ddmgr_->bddVar(Xs_[i]->indBddVars_[j][0]);
                    ind++;
                }
            }
            BDD* cubeCoarser = new BDD;
            *cubeCoarser = ddmgr_->bddComputeCube(varsCoarser, NULL, *system_->dimX_ - ones);
            cubesCoarser_.push_back(cubeCoarser);
            delete[] varsCoarser;
        }

        for (int i = 1; i < *system_->numAbs_; i++) {
            int* permuteCoarser = new int[numBDDVars_];
            for (int ind = 0; ind < numBDDVars_; ind++) {
                permuteCoarser[ind] = 0;
            }

            for (int dim = 0; dim < *system_->dimX_; dim++) {
                if (system_->etaRatio_[dim] == 2) {
                    for (size_t projInd = 1; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteCoarser[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd-1];
                    }
                }
                else if (system_->etaRatio_[dim] == 1){
                    for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteCoarser[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd];
                    }
                }
            }
            permutesCoarser_.push_back(permuteCoarser);

            clog << "permuteCoarser " << i << " to " << i-1 << ": ";
            printArray(permuteCoarser, numBDDVars_);
        }

        for (int i = 0; i < *system_->numAbs_ - 1; i++) {
            int* permuteFiner = new int[numBDDVars_];
            for (int ind = 0; ind < numBDDVars_; ind++) {
                permuteFiner[ind] = 0;
            }

            for (int dim = 0; dim < *system_->dimX_; dim++) {
                if (system_->etaRatio_[dim] == 2) {
                    for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteFiner[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i+1]->indBddVars_[dim][projInd+1];
                    }
                }
                else if (system_->etaRatio_[dim] == 1) {
                    for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteFiner[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i+1]->indBddVars_[dim][projInd];
                    }
                }
            }
            permutesFiner_.push_back(permuteFiner);

            clog << "permuteFiner " << i << " to " << i+1 << ": ";
            printArray(permuteFiner, numBDDVars_);
        }
    }

    /*! Initializes the BDDs that are the precursors to notXvars_ and notXUvars_ */
    void initializeCubes() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            BDD* cubeX = new BDD;
            BDD* cubeX2 = new BDD;

            BDD* varsX = new BDD[Xs_[i]->nvars_];
            BDD* varsX2 = new BDD[X2s_[i]->nvars_];

            for (size_t j = 0; j < Xs_[i]->nvars_; j++) {
                varsX[j] = ddmgr_->bddVar(Xs_[i]->idBddVars_[j]);
            }

            for (size_t j = 0; j < X2s_[i]->nvars_; j++) {
                varsX2[j] = ddmgr_->bddVar(X2s_[i]->idBddVars_[j]);
            }

            *cubeX = ddmgr_->bddComputeCube(varsX, NULL, Xs_[i]->nvars_);
            *cubeX2 = ddmgr_->bddComputeCube(varsX2, NULL, X2s_[i]->nvars_);


            cubesX_.push_back(cubeX);
            cubesX2_.push_back(cubeX2);

            delete[] varsX;
            delete[] varsX2;
        }

        cubeU_ = new BDD;
        *cubeU_ = U_->getCube();

        clog << "Initialized cubes.\n";
    }

    /*! Initializes the arrays of BDD variable IDs that allow a BDD over an X domain to be projected to the identical BDD over the corresponding X2 domain.
     */
    void initializePermutes() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            int* permuteXtoX2 = new int[numBDDVars_];
            int* permuteX2toX = new int[numBDDVars_];
            for (int j = 0; j < numBDDVars_; j++) {
                permuteXtoX2[j] = 0;
                permuteX2toX[j] = 0;
            }

            for (size_t j = 0; j < Xs_[i]->nvars_; j++) {
                permuteXtoX2[Xs_[i]->idBddVars_[j]] = X2s_[i]->idBddVars_[j];
                permuteX2toX[X2s_[i]->idBddVars_[j]] = Xs_[i]->idBddVars_[j];
            }

            permutesXtoX2_.push_back(permuteXtoX2);
            permutesX2toX_.push_back(permuteX2toX);
        }

        clog << "Initialized permutes.\n";
    }
    /*! Saves and prints to log file some information related to the reachability/always-eventually specification. */
    void verifySave() {
        clog << "etaXs_:\n";
        for (size_t i = 0; i < etaXs_.size(); i++) {
            clog << "abstraction " << i << ": ";
            for (int j = 0; j < *system_->dimX_; j++) {
                clog << etaXs_[i][j] << " ";
            }
            clog << '\n';
        }
        clog << "tau_:\n";
        for (size_t i = 0; i < etaXs_.size(); i++) {
            clog << "abstraction " << i << ": " << *tau_[i] << '\n';
        }
        printVec(Xs_, "X");
        //printVec(Os_, "O");
        cout << "U:\n";
        U_->printInfo(1);
        //printVec(Gs_, "G");

        checkMakeDir("plotting");
        Xs_[0]->writeToFile("plotting/X.bdd");
        SafeInner_[*system_->numAbs_-1]->writeToFile("plotting/SafeInner.bdd");
        checkMakeDir("SafeSets");
        saveVec(SafeInner_,"SafeSets/S");
    }
    
    /*! Initialize the vector containing packets of the symbolicSet for parallel computation of the transition relation */
    void initializeParallel(int ab) {
        for (int iter = 0; iter < numPacket_; iter++) {
            SymbolicSet* D = new SymbolicSet(*Xs_[ab]);
            splitDs_.push_back(D);
        }
    }

    /*! Saves a snapshot of a controller and its domain into the sequence of final controllers and controller domains.
        \param[in] ab	0-index of the abstraction which the controller and controller domain that should be saved belong to.
    */
    void saveCZ(int ab) {
//        if (Zs_[ab]->symbolicSet_ == ddmgr_->bddZero()) {
//            return;
//        }

        SymbolicSet* C = new SymbolicSet(*Cs_[ab]);
        C->symbolicSet_ = Cs_[ab]->symbolicSet_;
        finalCs_.push_back(C);
        SymbolicSet* Z = new SymbolicSet(*Xs_[ab]);
//        Z->symbolicSet_ = Zs_[ab]->symbolicSet_;
        Z->symbolicSet_ = Cs_[ab]->symbolicSet_.ExistAbstract(*cubeU_);
        finalZs_.push_back(Z);
        finalAbs_.push_back(ab);
    }

	void saveCZ() {
		for (int ab = 0; ab <= *system_->numAbs_ - 1; ab++) {
			saveCZ(ab);
		}
	}
};
}

#endif /* AdaptAbsSafe_HH_ */
