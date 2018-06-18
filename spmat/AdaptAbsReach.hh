#ifndef ADAPTABSREACH_HH_
#define ADAPTABSREACH_HH_

#include <cstdio>
#include <vector>
#include <queue>
#include <memory>
#include <climits>
#include <utility>
#include <iostream>
#include <stdexcept>

#include "UniformGrid.hh"
#include "System.hh"
#include "Helper.hh"
#include "RungeKutta4.hh"
#include "TransitionFunction.hh"
#include "TicToc.hh"
#include "Abstraction.hh"
#include "WinningDomain.hh"
#include "StateTree.hh"
#include "InputOutput.hh"

using std::clog;
using std::freopen;
using std::string;
using std::vector;
using namespace helper;

namespace scots {

/** @cond **/
/* default parameters for the solve_reachability_game */
namespace params {
  auto avoid = [](const abs_type&, const UniformGrid ss) noexcept {return false;};
  static std::vector<double> value {};
}

// enum ReachResult {CONVERGEDVALID, CONVERGEDINVALID, NOTCONVERGED};

class AdaptAbsReach {
public:
    System* system_; /*!< Contains abstraction parameters. */
    vector<double*> etaXs_; /*!< *system_->numAbs_ x *system_->dimX_ matrix of state space grid spacings. */
    vector<double*> tau_; /*!< *system_->numAbs_ x 1 matrix of time steps. */
    vector<UniformGrid*> Xs_; /*!< The *system_->numAbs_ "pre" state space abstractions, coarsest (0) to finest. */
    StateTree* tree; /*!< The data structure that maps states across layers, and takes care of the projection of winning states. */

  std::function<bool(const int, const UniformGrid&)> avoid; /*!< Function which evaluates to true for unsafe (obstacle) states of a given layer. */
	// std::function<bool(const int, const UniformGrid&)> target; /*!< Function which evaluates to true for winning states of a given layer. */
  std::function<bool(const abs_type, const int)> target; /*!< Function which evaluates to true for winning states of a given layer. First argument: state, second argument: layer. */
  //vector<UniformGrid*> Zs_; /*!< Instance of *Xs_[i] containing winning states. */
    vector<UniformGrid*> X2s_; /*!< The *system_->numAbs_ "post" state space abstractions, coarsest (0) to finest. */
	UniformGrid* U_; /*!< The single input space abstraction. */
    vector<UniformGrid*> Cs_; /*!< Controller \subseteq *Xs_[i] x *U_. */
    //int numBDDVars_; /*!< Total number of BDD variables used by ddmgr_. */
    //vector<BDD*> cubesX_; /*!< *cubesX_[i] is a BDD with a single minterm of all 1s over the domain of *Xs_[i]. */
    //vector<BDD*> cubesX2_; /*!< *cubesX2_[i] is a BDD with a single minterm of all 1s over the domain of *X2s_[i]. */
    //BDD* cubeU_; /*!< A BDD with a single minterm of all 1s over the domain of the input SymbolicSet. */
    //vector<BDD*> notXUvars_; /*!< notXUvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and *U_ and 1 otherwise. */
    //vector<BDD*> notXvars_; /*!< notXvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and 1 otherwise. */
    //vector<BDD*> cubesCoarser_;
    vector<TransitionFunction*> Ts_; /*!< Ts_[i] stores the transition relation of abstraction i. */
    vector<TransitionFunction*> TTs_; /*!< TTs_[i] is Ts_[i] with X2s_[i] existentially abtracted. */
    //vector<int*> permutesXtoX2_; /*!< To transform a BDD over X variables to an equivalent one over X2 variables. */
    //vector<int*> permutesX2toX_; /*!< To transform a BDD over X2 variables to an equivalent one over X variables. */
    //vector<int*> permutesCoarser_;
    //vector<int*> permutesFiner_;
    vector<OdeSolver*> solvers_; /*!< ODE solvers (Runge-Katta approximation) for each abstraction time step. */

    vector<UniformGrid*> Gs_; /*!< Instance of *Xs_[i] containing goal states. */
    vector<UniformGrid*> validZs_; /*!< Contains winning states that act as savepoints. */
    vector<UniformGrid*> validCs_; /*!< Controllers that act as savepoints. */
    vector<UniformGrid*> finalCs_; /*!< Sequence of controllers that satisfy the specification. */
    vector<UniformGrid*> finalZs_; /*!< Sequence of domains of finalCs_. */
    vector<int> finalAbs_; /*!< Sequence of abstraction layer indices of finalCs_. */
    vector<UniformGrid*> Ds_; /*!< Instance of *Xs_[i] containing possible winning states. */
    vector<UniformGrid*> innerDs_;
    vector<UniformGrid*> computedDs_;
    vector<UniformGrid*> savedZs_;

    vector<TransitionFunction*> uTs_; /*!< Transition relations for exploring, i.e. all with coarsest space gridding but varying time sampling. */

    int minToGoCoarser_;
    int minToBeValid_;
    int earlyBreak_;
    int verbose_;

    int m_;
    int p_;

    double abTime_;
    double synTime_;
    int rec_depth;


    /*!	Constructor for an AdaptAbsReach object.
     *  \param[in]  logFile     Filename of program log.
     */
    AdaptAbsReach(char* logFile) {
        freopen(logFile, "w", stderr);
        clog << logFile << '\n';

        abTime_ = 0;
        synTime_ = 0;
        rec_depth = 0;
    }
    /*! Destructor for an AdaptAbsReach object. */
    ~AdaptAbsReach() {
        clog << "Abstraction construction time: " << abTime_ << " seconds.\n";
        clog << "Controller synthesis time: " << synTime_ << " seconds.\n";
        deleteVecArray(etaXs_);
        deleteVec(tau_);
        deleteVec(Xs_);
        //deleteVec(Os_);
        //deleteVec(Zs_);
		/*delete obstacle;
		delete winning;*/
        deleteVec(X2s_);
        delete U_;
        deleteVec(Cs_);
        /*deleteVec(cubesX_);
        deleteVec(cubesX2_);
        delete cubeU_;
        deleteVec(notXUvars_);
        deleteVec(notXvars_);
        deleteVec(cubesCoarser_);*/
        deleteVec(Ts_);
        deleteVec(TTs_);
        /*deleteVecArray(permutesXtoX2_);
        deleteVecArray(permutesX2toX_);
        deleteVecArray(permutesCoarser_);
        deleteVecArray(permutesFiner_);*/
        deleteVec(solvers_);
        deleteVec(Gs_);
        deleteVec(validZs_);
        deleteVec(validCs_);
        deleteVec(finalCs_);
        deleteVec(finalZs_);
        deleteVec(Ds_);
        deleteVec(innerDs_);
        deleteVec(computedDs_);
        deleteVec(savedZs_);
        deleteVec(uTs_);
        fclose(stderr);
        /*delete ddmgr_;*/
    }

    template<class sys_type, class rad_type, class X_type, class U_type>
    void onTheFlyReach(int p, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        m_ = p; // max. iterations for consecutive reachability for non-coarsest layers
        p_ = p; // coarsest layer uncontrollable-pred. parameter
        minToBeValid_ = 2;
        clog << "m: " << m_ << '\n';
        clog << "p: " << p_ << '\n';

        // start by synthesizing all uTs_
		UniformGrid* D = new UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[0]);
		Ds_[0] = D;
		//Ds_[0] = UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[0]);
		clog << "Ds_[0] initialized";
		//delete D;


        TicToc timer;
        timer.tic();
        computeExplorationAbstractions(sysNext, radNext, x, u);
        abTime_ += timer.toc();
		clog << "Exploration abstractions computed";

		checkMakeDir("uTs");
		saveVec(uTs_, "uTs/uTs");

        // Ts_[0] is uTs_[0] with obstacles removed
		cout << "Starting computation of transition relation of layer 0: \n";
		Abstraction<X_type, U_type> abs(*Ds_[0], *U_); // coarsest state gridding
		abs.compute_gb(*Ts_[0], sysNext, radNext, *solvers_[0], avoid);

    // Non-lazy trial version : pre-compute all transitions
    for (int ab=1; ab < *system_->numAbs_; ab++) {
      UniformGrid* D = new UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[ab]);
  		Ds_[ab] = D;
      Abstraction<X_type, U_type> abs(*Ds_[ab], *U_); // coarsest state gridding
  		abs.compute_gb(*Ts_[ab], sysNext, radNext, *solvers_[ab], avoid);
    }
    // end
    // debug: testing reach
    // WinningDomain w;
    // w = reach(0, m_);
    // ReachResult result = w.get_result();
    // if (result == CONVERGEDVALID) {
    //   cout << "result: converged valid\n";
    // }
    // else if (result == CONVERGEDINVALID) {
    //   cout << "result: converged invalid\n";
    // }
    // else {
    //   cout << "result: not converged\n";
    // }
    //  end


        // begin on-the-fly reachability synthesis
        int ab = 0;
        int print = 1;
        checkMakeDir("C");
        onTheFlyReachRecurse(ab, sysNext, radNext, x, u, print);

        //clog << "controllers: " << finalCs_.size() << '\n';

        //checkMakeDir("C");
        //saveVec(finalCs_, "C/C");
        //checkMakeDir("Z");
        //saveVec(finalZs_, "Z/Z");
        //checkMakeDir("T");
        //saveVec(Ts_, "T/T");
        //clog << "Wrote Ts_ to file.\n";
        return;
    }

    template<class sys_type, class rad_type, class X_type, class U_type>
    void computeExplorationAbstractions(sys_type sysNext, rad_type radNext, X_type x, U_type u) {
		cout << "Starting computation of auxiliary abstractions:\n";
        for (int ab = 0; ab < *system_->numAbs_; ab++) {
			/*TransitionFunction	uTs_[ab];*/
			//// debug purpose
			//TransitionFunction uT; // exploration transition relations are all with coarsest gridding
			//for (int i = 0; i < *system_->numAbs_; i++) {
			//	uTs_.push_back(&uT);
			//}
			//// debug purpose
			cout << "Layer " << ab << ":\n";
            Abstraction<X_type, U_type> abs(*Ds_[0], *U_); // coarsest state gridding
            abs.compute_gb(*uTs_[ab], sysNext, radNext, *solvers_[ab]); // use solver with time step corresponding to that of each layer
        }
		cout << "Finished computation of auxiliary abstractions.\n\n";
    }

   template<class sys_type, class rad_type, class X_type, class U_type>
   void onTheFlyReachRecurse(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u, int print = 0) {
       rec_depth += 1;
       clog << '\n';
       clog << "current recursion depth: " << rec_depth << '\n';
       clog << "current abstraction: " << ab << '\n';
       clog << "controllers: " << finalCs_.size() << '\n';

       if (print) {
           cout << '\n';
           cout << "current abstraction: " << ab << '\n';
           cout << "controllers: " << finalCs_.size() << '\n';
       }

       WinningDomain w;
       if (ab == 0) {
           TicToc timer;
           timer.tic();
           w = reach(ab);
           synTime_ += timer.toc();
       }
       else {
           TicToc timer;
           timer.tic();
           w = reach(ab, m_);
           synTime_ += timer.toc();
       }
       ReachResult result = w.get_result();
       if (result == CONVERGEDVALID) {
           clog << "result: converged valid\n";
           if (print)
               cout << "result: converged valid\n";
       }
       else if (result == CONVERGEDINVALID) {
           clog << "result: converged invalid\n";
           if (print)
               cout << "result: converged invalid\n";
       }
       else if (result == NOTCONVERGED) {
           clog << "result: not converged\n";
           if (print)
               cout << "result: not converged\n";
       }

       if (result != CONVERGEDINVALID) {
           std::string file = "C/C" + std::to_string(rec_depth);
           write_to_file(StaticController(*Xs_[ab], *U_, std::move(w)), file);
           // validZs_[ab]->symbolicSet_ = Zs_[ab]->symbolicSet_;
           // validCs_[ab]->symbolicSet_ = Cs_[ab]->symbolicSet_;
           std::vector<abs_type> dom = w.get_winning_domain();
           for (abs_type i = 0; i < dom.size(); i++) {
             tree->markNode(ab, dom[i]);
           }
           clog << "saved as snapshot\n";
           if (print)
               cout << "saved as snapshot\n";
       }
       else { // result is CONVERGEDINVALID
           // Zs_[ab]->symbolicSet_ = validZs_[ab]->symbolicSet_; // reset this layer's progress
           // Cs_[ab]->symbolicSet_ = validCs_[ab]->symbolicSet_;
           clog << "reset to valids\n";
           if (print)
               cout << "reset to valids\n";
       }
       //
       if (result != NOTCONVERGED) { // ab = 0 always converges
           if (ab == *system_->numAbs_ - 1) {
               return;
           }
           else { // go finer
               clog << "Going finer\n";
               if (print)
                   cout << "Going finer\n";
               int nextAb = ab + 1;
               TicToc timer;
               // timer.tic();
               // eightToTen(ab, nextAb, sysNext, radNext, x, u);
               // abTime_ += timer.toc();
               // finer(Zs_[ab], Zs_[nextAb], ab);
               // Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
               // validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
               onTheFlyReachRecurse(nextAb, sysNext, radNext, x, u);
               return;
           }
       }
       else { // not converged, go coarser
           clog << "Going coarser\n";
           if (print)
               cout << "Going coarser\n";
           int nextAb = ab - 1;
           // if (nextAb != 0) { // pointless to do for coarsest layer
           //     TicToc timer;
           //     timer.tic();
           //     eightToTen(ab, nextAb, sysNext, radNext, x, u);
           //     abTime_ += timer.toc();
           // }
           // coarserInner(Zs_[nextAb], Zs_[ab], nextAb);
           // Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
           // validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
           onTheFlyReachRecurse(nextAb, sysNext, radNext, x, u);
           return;
       }
   }
//
//    template<class sys_type, class rad_type, class X_type, class U_type>
//    void eightToTen(int curAb, int nextAb, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
//        Ds_[curAb]->symbolicSet_ = validZs_[curAb]->symbolicSet_ | Gs_[curAb]->symbolicSet_; // target for uncontrollable reach
//        innerDs_[curAb]->symbolicSet_ = validZs_[curAb]->symbolicSet_ | Gs_[curAb]->symbolicSet_;
//        for (int c = curAb - 1; c >= 0; c--) {
//            coarserOuter(Ds_[c], Ds_[c+1], c); // compounding outer approximations?
//            coarserInner(innerDs_[c], innerDs_[c+1], c);
//        }
//        Ds_[0]->symbolicSet_ = uReach(Ds_[0]->symbolicSet_, curAb, p_); // do cooperative reach
//        Ds_[0]->symbolicSet_ &= (!innerDs_[0]->symbolicSet_); // states to do abstraction for
//        for (int c = 0; c < nextAb; c++) {
//            finer(Ds_[c], Ds_[c+1], c);
//        }
//        computeAbstraction(nextAb, sysNext, radNext, x, u);
//    }
//

    /**
     * fucntion: reach
     * @brief solve reachability game according to Algorithm 2 in  <a href="./../../manual/manual.pdf">manual</a>
     *
     * @param[in] trans_function - TransitionFunction of the symbolic model
     * @param[in] target - lambda expression of the form
     *                      \verbatim [] (const abs_type &i) -> bool \endverbatim
     *                      returns true if state i is in target set and false otherwise
     *
     * @param[in] avoid  - OPTIONALLY provide lambda expression of the form
     *                      \verbatim [] (const abs_type &i, const UniformGrid ss) -> bool \endverbatim
     *                      returns true if state i is in avoid set and false otherwise
     *
     * @param[in] m - OPTIONALLY provide maximum m of search for the reachability fix point.
     *                    Use m = -1 for the normal algorithm which runs until convergence.
     *                    The default value is -1
     *
     * @param[out] value - OPTIONALLY provide std::vector<double> value to obtain the value function
     *
     * @return -  WinningDomain that contains the set of winning states and valid inputs
     **/
    // template<class F1, class F2=decltype(params::avoid)>
    WinningDomain reach(//const TransitionFunction& trans_function,
                                          const int ab,
                                          // F1& target,
                                          // F2& avoid = params::avoid,
                                          const double m = -1,
                                          std::vector<double> & value = params::value ) {
      /* transition function of layer ab */
      TransitionFunction& trans_function = *Ts_[ab];
      /* state space of layer ab */
      UniformGrid ss = *Xs_[ab];
      /* target and obstacle regions */
      // F1& target = winning;
      // F2& avoid = obstacle;
      /* size of state alphabet */
      abs_type N=trans_function.m_no_states;
      /* size of input alphabet */
      abs_type M=trans_function.m_no_inputs;

      /* used to encode that a state is not in the winning domain */
      abs_type loosing = std::numeric_limits<abs_type>::max();
      if(M > loosing-1) {
        throw std::runtime_error("scots::reach: Number of inputs exceeds maximum supported value");
      }
      /* win_domain[i] = j
       * contains the input j associated with state i
       *
       * j = loosing if the target is not reachable from i
       *
       * initialize all states to loosing (win_domain[i]=loosing) */
      std::vector<abs_type> win_domain(N,loosing);
      /* initialize value */
      value.resize(N,std::numeric_limits<double>::infinity());
      /* keep track of the number of processed post */
      std::unique_ptr<abs_type[]> K(new abs_type[N*M]);
      /* keep track of the values (corresponds to M in Alg.2)*/
      std::unique_ptr<double[]>  edge_val(new double[N*M]);

      /* init fifo */
      std::queue<abs_type> fifo;
      for(abs_type i=0; i<N; i++) {
        // debug code begin
        std::cout << "state = " << i << ": ";
        std::array<double, 3> x;
        ss.itox(i, x);
        std::cout << "{ " << x[0] << ", " << x[1] << ", " << x[2] << " }" << " = ";
        if(target(i,ab))
          std::cout << "1, ";
        else
          std::cout << "0, ";

        if(avoid(i,ss))
          std::cout << "1" << '\n';
        else
          std::cout << "0" << '\n';
        // debug code end

        if(target(i, ab) && !avoid(i, ss)) {
          win_domain[i]=loosing;
          /* value is zero */
          value[i]=0;
          /* states in the target are added to the fifo */
          fifo.push(i);
        }
        for(abs_type j=0; j<M; j++) {
          edge_val[i*M+j]=0;
          K[i*M+j]=trans_function.m_no_post[i*M+j];
        }
      }

      ReachResult result = CONVERGEDINVALID;
      double lastValue;
      /* main loop */
      while(!fifo.empty()) {
        /* get state to be processed */
        abs_type q=fifo.front();
        fifo.pop();
        /* If the value of the popped state is greater or equal to m, break */
        if (value[q] > m && m != -1) {
          result = NOTCONVERGED;
          break;
        }
        /* Save the value of the popped state for comparison with minToBeValid_ */
        lastValue = value[q];
        /* loop over each input */
        for(abs_type j=0; j<M; j++) {
          /* loop over pre's associated with this input */
          for(abs_ptr_type v=0; v<trans_function.m_no_pre[q*M+j]; v++) {
            abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[q*M+j]+v];
            if(avoid(i, ss))
              continue;
            /* (i,j,q) is a transition */
            /* update the number of processed posts */
            K[i*M+j]--;
            /* update the max value of processed posts */
            edge_val[i*M+j]=(edge_val[i*M+j]>=1+value[q] ? edge_val[i*M+j] : 1+value[q]);
            /* check if for node i and input j all posts are processed */
            if(!K[i*M+j] && value[i]>edge_val[i*M+j]) {
              fifo.push(i);
              value[i]=edge_val[i*M+j];
              win_domain[i]=j;
            }
          }  /* end loop over all pres of state i under input j */
        }  /* end loop over all input j */
      }  /* fifo is empty */
      if (result != NOTCONVERGED && lastValue >= minToBeValid_) {
        /* if fifo is empty and the m reached is greater than minToBeValid_ */
        result = CONVERGEDVALID;
      }
      /* if the default value function was used, free the memory of the static object*/
      if(value == scots::params::value){
          value.clear();
          value.shrink_to_fit();
      }

      return WinningDomain(N,M,std::move(win_domain),std::vector<bool>{},result,loosing);
    }



//
//    // testing computeAbstraction
//    void test() {
//        Ds_[0]->addGridPoints();
//    }
//
//    // abstraction synthesis will only iterate over the elements in the domain of the state space SymbolicSet (first argument in constructing abstraction)
//    template<class sys_type, class rad_type, class X_type, class U_type>
//    void computeAbstraction(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
//        BDD D = Ds_[ab]->symbolicSet_ & (!computedDs_[ab]->symbolicSet_); // don't repeat sampling
//        D &= !Os_[ab]->symbolicSet_;
//        computedDs_[ab]->symbolicSet_ |= Ds_[ab]->symbolicSet_; // update computed part of transition relation
//        Ds_[ab]->symbolicSet_ = D;
//
//        if (ab == 1) {
//            checkMakeDir("D1");
//            if (!(std::ifstream("D1/1.bdd"))) {
//                Ds_[1]->writeToFile("D1/1.bdd");
//            } else if (!(std::ifstream("D1/2.bdd"))) {
//                Ds_[1]->writeToFile("D1/2.bdd");
//            } else if (!(std::ifstream("D1/3.bdd"))) {
//                Ds_[1]->writeToFile("D1/3.bdd");
//            } else if (!(std::ifstream("D1/4.bdd"))) {
//                Ds_[1]->writeToFile("D1/4.bdd");
//            } else if (!(std::ifstream("D1/5.bdd"))) {
//                Ds_[1]->writeToFile("D1/5.bdd");
//            } else if (!(std::ifstream("D1/6.bdd"))) {
//                Ds_[1]->writeToFile("D1/6.bdd");
//            }
//        }
//        if (ab == 2) {
//            checkMakeDir("D2");
//            if (!(std::ifstream("D2/1.bdd"))) {
//                Ds_[2]->writeToFile("D2/1.bdd");
//            } else if (!(std::ifstream("D2/2.bdd"))) {
//                Ds_[2]->writeToFile("D2/2.bdd");
//            } else if (!(std::ifstream("D2/3.bdd"))) {
//                Ds_[2]->writeToFile("D2/3.bdd");
//            } else if (!(std::ifstream("D2/4.bdd"))) {
//                Ds_[2]->writeToFile("D2/4.bdd");
//            } else if (!(std::ifstream("D2/5.bdd"))) {
//                Ds_[2]->writeToFile("D2/5.bdd");
//            } else if (!(std::ifstream("D2/6.bdd"))) {
//                Ds_[2]->writeToFile("D2/6.bdd");
//            }
//        }
//
//        SymbolicModelGrowthBound<X_type, U_type> abstraction(Ds_[ab], U_, X2s_[ab]);
//        abstraction.computeTransitionRelation(sysNext, radNext, *solvers_[ab]); // was hard-coded to 0, source of "tunneling" bug
//        if (abstraction.transitionRelation_ != ddmgr_->bddZero()) { // no point adding/displaying if nothing was added
//            Ts_[ab]->symbolicSet_ |= abstraction.transitionRelation_; // add to transition relation
////            BDD O2 = Os_[ab]->symbolicSet_.Permute(permutesXtoX2_[ab]); // causes unsound behavior, leaving here as warning
////            Ts_[ab]->symbolicSet_ &= !O2;
////            Ts_[ab]->printInfo(1);
//            TTs_[ab]->symbolicSet_ = Ts_[ab]->symbolicSet_.ExistAbstract(*notXUvars_[ab]);
//        }
//    }

    /*! Calculates all possible winning states at most p transitions away from Z.
     */
    /*BDD uReach(BDD Z, int ab, int p) {
        int i = 1;
        while (1) {
            BDD uPreZ = uPre(Z, ab);
            BDD N = uPreZ & (!Z);
            Z = uPreZ;

            if (i == p || N == ddmgr_->bddZero()) {
                return Z;
            }
            i += 1;
        }
    }*/

    /*! Calculates the uncontrollable predecessor of a given set with respect to the un-restricted transition relation at the coarsest level of abstraction.
     *  \param[in]  Z       The given set.
     *  \return     BDD containing {x} for which there exists a u s.t. there exists a post state of (x,u) in Z.
     */
    /*BDD uPre(BDD Z, int ab) {
        BDD Z2 = Z.Permute(permutesXtoX2_[0]);
        BDD uPreZ = uTs_[ab]->symbolicSet_.AndAbstract(Z2, *notXUvars_[0]);
        uPreZ = uPreZ.ExistAbstract(*notXvars_[0]);
        return uPreZ;
    }*/


    /*! Calculates the enforceable predecessor of the given set with respect to the transition relation at the specified level of abstraction.
     *  \param[in]  Z           The winning set.
     *  \param[in]  curAbs      0-index of the current abstraction.
     *  \return     BDD containing {(x,u)} for which all post states are in Z.
     */
    //BDD cPre(BDD Z, int curAbs) {
    //    // swap to X2s_[curAbs]
    //    BDD Z2 = Z.Permute(permutesXtoX2_[curAbs]);
    //    // posts outside of winning set
    //    BDD nZ2 = !Z2;
    //    // {(x,u)} with some posts outside of winning set
    //    BDD Fbdd = Ts_[curAbs]->symbolicSet_.AndAbstract(nZ2, *notXUvars_[curAbs]);
    //    // {(x,u)} with no posts outside of winning set
    //    BDD nF = !Fbdd;
    //    // get rid of junk
    //    BDD preF = TTs_[curAbs]->symbolicSet_.AndAbstract(nF, *notXUvars_[curAbs]);
    //    return preF;
    //}

    /*void finer(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
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
    }*/

    /*! Initializes data members.
     *  \param[in]  system      Contains abstraction parameters.
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     */
    template<class isG, class isO>
    void initialize(System* system, isG G, isO O) {
        //ddmgr_ = new Cudd;
        system_ = system;

        TicToc timer;
        timer.tic();

        initializeEtaTau();
        initializeSolvers();
        initializeSymbolicSets();
        initializeFunctions(G, O);
        /*initializeNumBDDVars();
        initializePermutes();
        initializeCubes();
        initializeNotVars();
        initializeProjs();*/
        //verifySave();

        abTime_ += timer.toc();
    }
    template<class isG, class isO>
    void initializeFunctions(isG G, isO O) {
      // target = G;
      avoid = O;
      int numAbs = *system_->numAbs_;
      abs_type N = Xs_[numAbs-1]->size();
      for (abs_type i = 0; i < N; i++) {
        if (G(i, Xs_[numAbs-1])) {
          tree->markNode(numAbs-1, i);
        }
      }
      auto target = [&](const abs_type &abs_state, const int ab) {
        if (tree->getMarkingStatus(ab, abs_state))
          return true;
        else
          return false;
      };
    }

    /*! Initializes SymbolicSet data members.
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     *  \param[in]	addG	Function pointer specifying the points that should be added to the goal set.
     */
    void initializeSymbolicSets() {
        for (int i = 0; i < *system_->numAbs_; i++) {
			UniformGrid* X = new UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[i]);
			Xs_.push_back(X);
			//X2s_.push_back(X);
        }
		clog << "Xs_ and X2s_ initialized with full domain.\n";

      tree = new StateTree(*system_->numAbs_, Xs_, system_->etaRatio_);

        U_ = new UniformGrid(*system_->dimU_, system_->lbU_, system_->ubU_, system_->etaU_);
        clog << "U_ initialized with full domain.\n";
		//delete U;

		/*O_type obstacle = isO;
		G_type winning = isG;*/
        //clog << "Os_ initialized according to specification.\n";
		//clog << "Obstacle and winning functions initialized according to the specification.\n";

        /*for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* Z = new SymbolicSet(*Xs_[i]);
            Zs_.push_back(Z);

        }
        clog << "Zs_ initialized with empty domain.\n";*/

        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    SymbolicSet* C = new SymbolicSet(*Xs_[i], *U_);
        //    Cs_.push_back(C);
        //}
        //clog << "Cs_ initialized with empty domain.\n";

        //initializeSpec(&Gs_, addG);
        //clog << "Gs_ initialized according to specification.\n";

        //// checking that specification is valid
        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    if ((Gs_[i]->symbolicSet_ & Os_[i]->symbolicSet_) != ddmgr_->bddZero()) {
        //        error("Error: G and O have nonzero intersection.\n");
        //    }
        //}
        //clog << "No obstacle problem with specification.\n";

        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    SymbolicSet* validZ = new SymbolicSet(*Xs_[i]);
        //    validZs_.push_back(validZ);
        //}
        //clog << "validZs_ initialized with empty domain.\n";

        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    SymbolicSet* validC = new SymbolicSet(*Xs_[i], *U_);
        //    validCs_.push_back(validC);
        //}
        //clog << "validCs_ initialized with empty domain.\n";

		for (int i = 0; i < *system_->numAbs_; i++) {
			UniformGrid* D = new UniformGrid;
            Ds_.push_back(D);
        }
        clog << "Ds_ initialized with empty domain.\n";
		//delete D;

        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    SymbolicSet* innerD = new SymbolicSet(*Xs_[i]);
        //    innerDs_.push_back(innerD);
        //}
        //clog << "innerDs_ initialized with empty domain.\n";

        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    SymbolicSet* computedD = new SymbolicSet(*Xs_[i]);
        //    computedDs_.push_back(computedD);
        //}
        //clog << "computedDs_ initialized with empty domain.\n";

        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    SymbolicSet* savedZ = new SymbolicSet(*Xs_[i]);
        //    savedZs_.push_back(savedZ);
        //}
        //clog << "savedZs_ initialized with empty domain.\n";

		for (int i = 0; i < *system_->numAbs_; i++) {
			TransitionFunction* Ts = new TransitionFunction; // exploration transition relations are all with coarsest gridding
			Ts_.push_back(Ts);
		}
		clog << "Ts_ initialized with empty domain.\n";

        //for (int i = 0; i < *system_->numAbs_; i++) {
        //    SymbolicSet* TT = new SymbolicSet(*Cs_[i]);
        //    TTs_.push_back(TT);
        //}
        //clog << "TTs_ initialized with empty domain.\n";

		for (int i = 0; i < *system_->numAbs_; i++) {
			TransitionFunction* uT = new TransitionFunction; // exploration transition relations are all with coarsest gridding
            uTs_.push_back(uT);
        }
        clog << "uT_ initialized with empty domain.\n";
		//delete uT;
    }

    ///*! Initializes a vector of SymbolicSets that is an instance of Xs_ containing points as specified by addSpec.
    // *  \param[in] vec      Vector of SymbolicSets to initialize.
    // *  \param[in] addSpec  Function pointer specifying the points to add to each element of vec.
    // */
    //template<class vec_type, class spec_type>
    //void initializeSpec(vec_type* vec, spec_type addSpec) {
    //    for (int i = 0; i < *system_->numAbs_; i++) {
    //        SymbolicSet* Xinstance = new SymbolicSet(*Xs_[i]);
    //        addSpec(Xinstance);
    //        vec->push_back(Xinstance);
    //    }
    //}

    ///*! Initializes the BDDs useful for existential abstraction. Must be called only after initializing Xs, U, and X2s. */
    //void initializeNotVars() {
    //    for (int i = 0; i < *system_->numAbs_; i++) {
    //        BDD* notXUvars = new BDD;
    //        BDD* notXvars = new BDD;
    //        *notXUvars = ddmgr_->bddOne();
    //        for (int j = 0; j < *system_->numAbs_; j++) {
    //            *notXUvars &= *cubesX2_[j];
    //            if (i != j) {
    //                *notXUvars &= *cubesX_[j];
    //            }
    //        }
    //        *notXvars = *notXUvars & *cubeU_;

    //        notXUvars_.push_back(notXUvars);
    //        notXvars_.push_back(notXvars);
    //    }

    //    clog << "Initialized notVars.\n";
    //}

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

    ///*! Initializes the number of total distinct BDD variables in use. */
    //void initializeNumBDDVars() {
    //    numBDDVars_ = 0;
    //    for (int i = 0; i < *system_->numAbs_; i++) {
    //        numBDDVars_ += Xs_[i]->nvars_;
    //        numBDDVars_ += X2s_[i]->nvars_;
    //    }
    //    numBDDVars_ += U_->nvars_;

    //    clog << "Number of BDD variables: " << numBDDVars_ << '\n';
    //}

    //void initializeProjs() {
    //    int ones = 0;
    //    for (int i = 0; i < *system_->dimX_; i++) {
    //        if (system_->etaRatio_[i] == 1) {
    //            ones++;
    //        }
    //    }

    //    for (int i = 1; i < *system_->numAbs_; i++) {
    //        BDD* varsCoarser = new BDD[*system_->dimX_ - ones];
    //        int ind = 0;
    //        for (int j = 0; j < *system_->dimX_; j++) {
    //            if (system_->etaRatio_[j] == 2) {
    //                varsCoarser[ind] = ddmgr_->bddVar(Xs_[i]->indBddVars_[j][0]);
    //                ind++;
    //            }
    //        }
    //        BDD* cubeCoarser = new BDD;
    //        *cubeCoarser = ddmgr_->bddComputeCube(varsCoarser, NULL, *system_->dimX_ - ones);
    //        cubesCoarser_.push_back(cubeCoarser);
    //        delete[] varsCoarser;
    //    }

    //    for (int i = 1; i < *system_->numAbs_; i++) {
    //        int* permuteCoarser = new int[numBDDVars_];
    //        for (int ind = 0; ind < numBDDVars_; ind++) {
    //            permuteCoarser[ind] = 0;
    //        }

    //        for (int dim = 0; dim < *system_->dimX_; dim++) {
    //            if (system_->etaRatio_[dim] == 2) {
    //                for (size_t projInd = 1; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
    //                    permuteCoarser[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd-1];
    //                }
    //            }
    //            else if (system_->etaRatio_[dim] == 1){
    //                for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
    //                    permuteCoarser[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd];
    //                }
    //            }
    //        }
    //        permutesCoarser_.push_back(permuteCoarser);

    //        clog << "permuteCoarser " << i << " to " << i-1 << ": ";
    //        printArray(permuteCoarser, numBDDVars_);
    //    }

    //    for (int i = 0; i < *system_->numAbs_ - 1; i++) {
    //        int* permuteFiner = new int[numBDDVars_];
    //        for (int ind = 0; ind < numBDDVars_; ind++) {
    //            permuteFiner[ind] = 0;
    //        }

    //        for (int dim = 0; dim < *system_->dimX_; dim++) {
    //            if (system_->etaRatio_[dim] == 2) {
    //                for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
    //                    permuteFiner[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i+1]->indBddVars_[dim][projInd+1];
    //                }
    //            }
    //            else if (system_->etaRatio_[dim] == 1) {
    //                for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
    //                    permuteFiner[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i+1]->indBddVars_[dim][projInd];
    //                }
    //            }
    //        }
    //        permutesFiner_.push_back(permuteFiner);

    //        clog << "permuteFiner " << i << " to " << i+1 << ": ";
    //        printArray(permuteFiner, numBDDVars_);
    //    }
    //}

    ///*! Initializes the BDDs that are the precursors to notXvars_ and notXUvars_ */
    //void initializeCubes() {
    //    for (int i = 0; i < *system_->numAbs_; i++) {
    //        BDD* cubeX = new BDD;
    //        BDD* cubeX2 = new BDD;

    //        BDD* varsX = new BDD[Xs_[i]->nvars_];
    //        BDD* varsX2 = new BDD[X2s_[i]->nvars_];

    //        for (size_t j = 0; j < Xs_[i]->nvars_; j++) {
    //            varsX[j] = ddmgr_->bddVar(Xs_[i]->idBddVars_[j]);
    //        }

    //        for (size_t j = 0; j < X2s_[i]->nvars_; j++) {
    //            varsX2[j] = ddmgr_->bddVar(X2s_[i]->idBddVars_[j]);
    //        }

    //        *cubeX = ddmgr_->bddComputeCube(varsX, NULL, Xs_[i]->nvars_);
    //        *cubeX2 = ddmgr_->bddComputeCube(varsX2, NULL, X2s_[i]->nvars_);


    //        cubesX_.push_back(cubeX);
    //        cubesX2_.push_back(cubeX2);

    //        delete[] varsX;
    //        delete[] varsX2;
    //    }

    //    cubeU_ = new BDD;
    //    *cubeU_ = U_->getCube();

    //    clog << "Initialized cubes.\n";
    //}

    ///*! Initializes the arrays of BDD variable IDs that allow a BDD over an X domain to be projected to the identical BDD over the corresponding X2 domain.
    // */
    //void initializePermutes() {
    //    for (int i = 0; i < *system_->numAbs_; i++) {
    //        int* permuteXtoX2 = new int[numBDDVars_];
    //        int* permuteX2toX = new int[numBDDVars_];
    //        for (int j = 0; j < numBDDVars_; j++) {
    //            permuteXtoX2[j] = 0;
    //            permuteX2toX[j] = 0;
    //        }

    //        for (size_t j = 0; j < Xs_[i]->nvars_; j++) {
    //            permuteXtoX2[Xs_[i]->idBddVars_[j]] = X2s_[i]->idBddVars_[j];
    //            permuteX2toX[X2s_[i]->idBddVars_[j]] = Xs_[i]->idBddVars_[j];
    //        }

    //        permutesXtoX2_.push_back(permuteXtoX2);
    //        permutesX2toX_.push_back(permuteX2toX);
    //    }

    //    clog << "Initialized permutes.\n";
    //}
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
        //printVec(Xs_, "X");
        /*printVec(Os_, "O");*/
        cout << "U:\n";
        //U_->printInfo(1);
        /*printVec(Gs_, "G");*/

        checkMakeDir("plotting");
        //Xs_[0]->writeToFile("plotting/X.bdd");
        /*Os_[*system_->numAbs_-1]->writeToFile("plotting/O.bdd");
        Gs_[*system_->numAbs_-1]->writeToFile("plotting/G.bdd");
        checkMakeDir("O");
        saveVec(Os_, "O/O");
        checkMakeDir("G");
        saveVec(Gs_, "G/G");*/
    }

    ///*! Saves a snapshot of a controller and its domain into the sequence of final controllers and controller domains.
    //    \param[in] ab	0-index of the abstraction which the controller and controller domain that should be saved belong to.
    //*/
    //void saveCZ(int ab) {
    //    if (Zs_[ab]->symbolicSet_ == ddmgr_->bddZero()) {
    //        return;
    //    }

    //    SymbolicSet* C = new SymbolicSet(*Cs_[ab]);
    //    C->symbolicSet_ = Cs_[ab]->symbolicSet_;
    //    finalCs_.push_back(C);
    //    SymbolicSet* Z = new SymbolicSet(*Xs_[ab]);
    //    Z->symbolicSet_ = Zs_[ab]->symbolicSet_;
    //    finalZs_.push_back(Z);
    //    finalAbs_.push_back(ab);
    //}
};
}

#endif /* ADAPTABSREACH_HH_ */
