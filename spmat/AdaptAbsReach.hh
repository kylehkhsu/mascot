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
#include "Helper_BDD.hh"
#include "RungeKutta4.hh"
#include "TransitionFunction.hh"
#include "TicToc.hh"
#include "Abstraction.hh"
#include "WinningDomain.hh"
#include "StateTree.hh"
#include "InputOutput.hh"
#include "Goal.hh"

using std::clog;
using std::freopen;
using std::string;
using std::vector;
using namespace helper;

namespace scots {

/** @cond **/
/* default parameters for the solve_reachability_game */
namespace params {
  auto avoid = [](const abs_type&, const UniformGrid* ss) noexcept {return false;};
  static std::vector<double> value {};
}

// enum ReachResult {CONVERGEDVALID, CONVERGEDINVALID, NOTCONVERGED};

class AdaptAbsReach {
public:
    System* system_; /*!< Contains abstraction parameters. */
    vector<double*> etaXs_; /*!< *system_->numAbs_ x *system_->dimX_ matrix of state space grid spacings. */
    vector<double*> tau_; /*!< *system_->numAbs_ x 1 matrix of time steps. */
    vector<UniformGrid*> Xs_; /*!< The *system_->numAbs_ "pre" state space abstractions, coarsest (0) to finest. */
    StateTree* tree_; /*!< The data structure that maps states across layers, and takes care of the projection of winning states. */

  std::function<bool(const int, const UniformGrid*)> avoid; /*!< Function which evaluates to true for unsafe (obstacle) states of a given layer. */
	// std::function<bool(const int, const UniformGrid&)> target; /*!< Function which evaluates to true for winning states of a given layer. */
  // std::function<bool(const abs_type, const int)> target; /*!< Function which evaluates to true for winning states of a given layer. First argument: state, second argument: layer. */
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

    vector<Goal*> Gs_; /*!< Instance of *Xs_[i] containing goal states. */
    vector<UniformGrid*> validZs_; /*!< Contains winning states that act as savepoints. */
    vector<UniformGrid*> validCs_; /*!< Controllers that act as savepoints. */
    vector<StaticController*> finalCs_; /*!< Sequence of controllers that satisfy the specification. */
    // vector<vector<abs_type>*> finalZs_; /*!< Sequence of domains of finalCs_. */
    vector<Goal*> finalZs_; /*!< Sequence of domains of finalCs_. */
    vector<int> finalAbs_; /*!< Sequence of abstraction layer indices of finalCs_. */
    vector<UniformGrid*> Ds_; /*!< Instance of *Xs_[i] containing possible winning states. */
    vector<UniformGrid*> innerDs_;
    vector<UniformGrid*> computedDs_;
    vector<UniformGrid*> savedZs_;

    // debug purpose: printing the intermediate states_to_explore
    std::vector<Goal*> intermediate_states_to_explore;
    // end

    vector<TransitionFunction*> uTs_; /*!< Transition relations for exploring, i.e. all with coarsest space gridding but varying time sampling. */

    int minToGoCoarser_;
    int minToBeValid_;
    int earlyBreak_;
    bool verbose_;

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
        // for (size_t i = 0; i < Gs_.size(); i++) {
        //   Gs_[i]->clear();
        // }
        // Gs_.clear();
        deleteVec(validZs_);
        deleteVec(validCs_);
        deleteVec(finalCs_);
        deleteVec(finalZs_);
        // for (size_t i = 0; i < finalZs_.size(); i++) {
        //   finalZs_[i]->clear();
        // }
        // finalZs_.clear();
        deleteVec(Ds_);
        deleteVec(innerDs_);
        deleteVec(computedDs_);
        deleteVec(savedZs_);
        deleteVec(uTs_);
        fclose(stderr);
        /*delete ddmgr_;*/
    }

    template<class sys_type, class rad_type, class X_type, class U_type>
    void onTheFlyReach(int p, sys_type sysNext, rad_type radNext, X_type x, U_type u, bool lazy=true, int readAbs=0) {
        m_ = p; // max. iterations for consecutive reachability for non-coarsest layers
        p_ = 2; // coarsest layer uncontrollable-pred. parameter
        minToBeValid_ = 2;
        clog << "m: " << m_ << '\n';
        clog << "p: " << p_ << '\n';

        UniformGrid* D = new UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[0]);
        Ds_[0] = D;
        //Ds_[0] = UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[0]);
        clog << "Ds_[0] initialized";
        //delete D;
        // Ts_[0] is uTs_[0] with obstacles removed
        checkMakeDir("T");
        Abstraction<X_type, U_type> abs(*Ds_[0], *U_); // coarsest state gridding
        if (readAbs==0) {
          cout << "Starting computation of transition relation of layer 0: \n";
          abs.compute_gb(*Ts_[0], sysNext, radNext, *solvers_[0], avoid);
          write_to_file(*Ts_[0], "T/T0", true);
          if (verbose_) {
            Ts_[0]->print_info();
            tree_->print_info();
          }
        }
        else {
          bool flag = read_from_file(*Ts_[0], "T/T0");
          if (!flag)
            throw std::runtime_error("\nAdaptAbsReach: could not read transition function from file");
          else {
            std::cout << "Read transition of layer 0 from file." << '\n';
            if (verbose_)
              Ts_[0]->print_info();
          }
        }

        /* debug purpose: printing the transition matrix */
        // Ts_[0]->print_domain(Xs_[0], U_);
        // const abs_ptr_type ntr = Ts_[0]->get_no_transitions();
        // double** tr_domain = new double*[ntr];
        // for(int i = 0; i < ntr; ++i)
        //     tr_domain[i] = new double[2*(*system_->dimX_)+ (*system_->dimU_)];
        //
        // Ts_[0]->get_domain(Xs_[0], U_, tr_domain);
        // printVecVec(tr_domain, "transition");
        //
        // delete[] tr_domain;

        // start by synthesizing all uTs_



        if (lazy) {
          TicToc timer;
          timer.tic();
          computeExplorationAbstractions(sysNext, radNext, x, u);
          abTime_ += timer.toc();
          clog << "Exploration abstractions computed";

          checkMakeDir("uTs");
          saveVec(uTs_, "uTs/uTs");
        }
        // Non-lazy trial version : pre-compute all transitions
        else {
          for (int ab=1; ab < *system_->numAbs_; ab++) {
            UniformGrid* D = new UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[ab]);
        		Ds_[ab] = D;
            Abstraction<X_type, U_type> abs(*Ds_[ab], *U_); // coarsest state gridding
            if (readAbs==0) {
              cout << "Starting computation of transition relation of layer " << ab << ": \n";
              abs.compute_gb(*Ts_[ab], sysNext, radNext, *solvers_[ab], avoid);
              std::string file = "T/T" + std::to_string(ab);
              write_to_file(*Ts_[ab], file);
              if (verbose_)
                Ts_[ab]->print_info();
            }
        	  else {
              std::string file = "T/T" + std::to_string(ab);
              bool flag = read_from_file(*Ts_[ab], file);
              if (!flag)
                throw std::runtime_error("\nAdaptAbsReach: could not read transition function from file");
              else {
                std::cout << "Read transition of layer " << ab << " from file." << '\n';
                if (verbose_)
                  Ts_[ab]->print_info();
              }
            }
          }
          /* debug purpose: use the game solver of scots v0.2 to synthesize controller for the finest layer */
          // WinningDomain win=solve_reachability_game(*system_->numAbs_-1);
          // std::cout << "Winning domain size: " << win.get_size() << std::endl;
          //
          // std::cout << "\nWrite controller to controller.scs \n";
          // if(write_to_file(scots::StaticController(*Xs_[*system_->numAbs_-1],*U_,std::move(win)),"controller"))
          // std::cout << "Done. \n";
          // return;
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
        onTheFlyReachRecurse(ab, sysNext, radNext, x, u, lazy, print);

        // debug purpose: reading controller from file
        // scots::StaticController *controller = new scots::StaticController;
        // if(!scots::read_from_file(*controller,"C/C1.scs")) {
        //   bool debug = false;
        // }

        clog << "controllers: " << finalCs_.size() << '\n';

        checkMakeDir("G");
        saveVec(Gs_, "G/G");
        checkMakeDir("C");
        saveVec(finalCs_, "C/C");
        checkMakeDir("Z");
        saveVec(finalZs_, "Z/Z");
        checkMakeDir("T");
        saveVec(Ts_, "T/T");
        // debug purpose: print intermediate_states_to_explore
        checkMakeDir("D");
        saveVec(intermediate_states_to_explore, "D/D");
        clog << "Wrote Ts_ to file.\n";

        /* debug purpose */
        // Goal* debug;
        // bool  flag = read_from_file(*debug, "Z/Z1");
        // StaticController* debug2;
        // bool flag = read_from_file(*debug2, "C/C1");

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
   void onTheFlyReachRecurse(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u, bool lazy = true, int print = 0) {
       rec_depth += 1;
       clog << '\n';
       clog << "current recursion depth: " << rec_depth << '\n';
       clog << "current abstraction: " << ab << '\n';
       // clog << "controllers: " << finalCs_.size() << '\n';

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

       // if ((result != CONVERGEDINVALID) && (ab!=*system_->numAbs_-1)) {
        // if ((result != CONVERGEDINVALID) || (ab==*system_->numAbs_-1)) {
       if (result != CONVERGEDINVALID) {
         // validZs_[ab]->symbolicSet_ = Zs_[ab]->symbolicSet_;
         // validCs_[ab]->symbolicSet_ = Cs_[ab]->symbolicSet_;
         std::vector<abs_type> dom = w.get_winning_domain();
         for (abs_type i = 0; i < dom.size(); i++) {
           tree_->markNode(ab, dom[i]);
         }
         saveCZ(ab, w);
         // std::cout << "Winning domain size: " << w.get_size() << std::endl;
         // std::string file = "C/C" + std::to_string(rec_depth);
         // write_to_file(StaticController(*Xs_[ab], *U_, std::move(w), *tau_[ab]), file);
         clog << "saved as snapshot\n";
         if (print)
             cout << "saved as snapshot\n";
       }
       else { // result is CONVERGEDINVALID
         // Zs_[ab]->symbolicSet_ = validZs_[ab]->symbolicSet_; // reset this layer's progress
         // Cs_[ab]->symbolicSet_ = validCs_[ab]->symbolicSet_;
         clog << "reset to valids\n"; // reset by not marking any new node in tree_
         if (print)
             cout << "reset to valids\n";
       }
       //
       if (result != NOTCONVERGED) { // ab = 0 always converges
           if (ab == *system_->numAbs_ - 1) {
             // std::vector<abs_type> dom = w.get_winning_domain();
             // saveCZ(ab, w);
             // for (abs_type i = 0; i < dom.size(); i++) {
             //   tree_->markNode(ab, dom[i]);
             // }
             // std::cout << "Winning domain size: " << w.get_size() << std::endl;
             // std::string file = "C/C" + std::to_string(rec_depth);
             // write_to_file(StaticController(*Xs_[ab], *U_, std::move(w)), file);
             return;
           }
           else { // go finer
               clog << "Going finer\n";
               if (print)
                   cout << "Going finer\n";
               int nextAb = ab + 1;
               // TicToc timer;
               // timer.tic();
               // eightToTen(ab, nextAb, sysNext, radNext, x, u);
               // abTime_ += timer.toc();
               // finer(Zs_[ab], Zs_[nextAb], ab);
               // Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
               // validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;

               if (lazy) {
                 TicToc timer;
                 timer.tic();
                 explore(nextAb, x, u, sysNext, radNext);
                 abTime_ += timer.toc();
               }
               onTheFlyReachRecurse(nextAb, sysNext, radNext, x, u, lazy, 1);
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

           if (lazy) {
             TicToc timer;
             timer.tic();
             explore(nextAb, x, u, sysNext, radNext);
             abTime_ += timer.toc();
           }
           onTheFlyReachRecurse(nextAb, sysNext, radNext, x, u, lazy, 1);
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
     * function: target
     * @brief check if a state is in the target or not
     *
     * @param[in] abs_state - the symbolic state
     * @param[in] ab - the abstraction layer
     **/
     bool target(const abs_type abs_state, const int ab) {
       if (tree_->getMarkingStatus(ab, abs_state)==2)
         return true;
       else
         return false;
     }

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
      abs_ptr_type no_trans = trans_function.get_no_transitions();
      /* state space of layer ab */
      UniformGrid* ss = Xs_[ab];
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
        // std::cout << "state = " << i << ": ";
        // std::array<double, 3> x;
        // ss.itox(i, x);
        // std::cout << "{ " << x[0] << ", " << x[1] << ", " << x[2] << " }" << " = ";
        // if(target(i,ab))
        //   std::cout << "1, ";
        // else
        //   std::cout << "0, ";
        //
        // if(avoid(i,ss))
        //   std::cout << "1" << '\n';
        // else
        //   std::cout << "0" << '\n';
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
      double lastValue=0;
      /* main loop */
      while(!fifo.empty()) {
        /* get state to be processed */
        abs_type q=fifo.front();
        fifo.pop();
        /* If the value of the popped state is greater or equal to m, break */
        if (value[q] >= m && m != -1) {
          result = NOTCONVERGED;
          break;
        }
        /* Save the value of the popped state for comparison with minToBeValid_ */
        lastValue = value[q];
        /* loop over each input */
        for(abs_type j=0; j<M; j++) {
          /* loop over pre's associated with this input */
          for(abs_ptr_type v=0; v<trans_function.m_no_pre[q*M+j]; v++) {
            std::vector<abs_type> debug = trans_function.get_pre(q, j);
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
      // debug purpose
      std::cout << "Fix point iterations: " << lastValue << '\n';
      /* if the default value function was used, free the memory of the static object*/
      if(value == scots::params::value){
          value.clear();
          value.shrink_to_fit();
      }

      return WinningDomain(N,M,std::move(win_domain),std::vector<bool>{},result,loosing);
    }

    template<class sys_type, class rad_type, class X_type, class U_type>
    void explore(int ab, X_type x, U_type u, sys_type sysNext, rad_type radNext) {
      /* for coarsest layer, do nothing */
      if (ab==0)
        return;
      /* auxiliary transition function of layer ab */
      TransitionFunction& aux_trans_function = *uTs_[ab];
      /* transition function of layer ab */
      TransitionFunction& trans_function = *Ts_[ab];
      /* state space of layer ab */
      UniformGrid ss = *Xs_[ab];
      /* size of input alphabet */
      abs_type M=aux_trans_function.m_no_inputs;
      int infinity = std::numeric_limits<int>::max() - 1;
      // int infinity = -1;
      /* vector containing the states to explore */
      std::vector<abs_type> states_to_explore;

      // int* val = NULL;
      // val = new int[Xs_[0]->size()];
      std::vector<int> val;
      for (abs_type i = 0; i < Xs_[0]->size(); i++) {
        val.push_back(infinity);
      }
      // std::vector<abs_type> states_to_explore;
      // abs_type states_to_explore[Xs_[0]->size()];
      for (abs_type q = 0; q < Xs_[0]->size(); q++) {
        if (tree_->getMarkingStatus(0, q) > 0) {
          val[q] = 0;
          // for(abs_type j=0; j<M; j++) {
          //   /* loop over pre's associated with this input */
          //   for(abs_ptr_type v=0; v<aux_trans_function.m_no_pre[q*M+j]; v++) {
          //     abs_type i=aux_trans_function.m_pre[aux_trans_function.m_pre_ptr[q*M+j]+v];
          //     // states_to_explore.push_back(i);
          //     val[i] = 1;
          //   }
          // }
        }
      }

      /* cooperative predecessor */
      for (int i = 0; i < p_; i++) {
        for (abs_type q = 0; q < Xs_[0]->size(); q++) {
          for (abs_type j=0; j<M; j++) {
            /* loop over pre's associated with this input */
            for(abs_ptr_type v=0; v<aux_trans_function.m_no_pre[q*M+j]; v++) {
              abs_type r=aux_trans_function.m_pre[aux_trans_function.m_pre_ptr[q*M+j]+v];
              // states_to_explore.push_back(i);
              // val[r] = (val[r] > val[q] + 1) ? (val[q] + 1) : (val[r] == -1 ? (val[q]+1) : val[r]);
              val[r] = (val[r] > val[q] + 1) ? (val[q] + 1) : val[r];
            }
          }
        }
      }


      // int index = 0;
      // while ((val[states_to_explore[index]] <= p_)
      //         || (index < states_to_explore.size())) {
      //   abs_type q = states_to_explore[index];
      //   for(abs_type j=0; j<M; j++) {
      //     /* loop over pre's associated with this input */
      //     for(abs_ptr_type v=0; v<aux_trans_function.m_no_pre[q*M+j]; v++) {
      //       abs_type i=aux_trans_function.m_pre[aux_trans_function.m_pre_ptr[q*M+j]+v];
      //       states_to_explore.push_back(i);
      //       val[i] = val[q] + 1;
      //     }
      //   }
      //   index += 1;
      // }

        /* projecting down */
      for (abs_type i = 0; i < Xs_[0]->size(); i++) {
        if (val[i] > 0 && val[i] <= p_)
          states_to_explore.push_back(i);
      }
      std::vector<abs_type> temp;
      for (int i = 0; i < ab; i++) {
        for (int j = 0; j < states_to_explore.size(); j++) {
          StateTreeNode* node = tree_->getNode(i, states_to_explore[j]);
          for (int k = 0; k < node->get_no_child(); k++) {
            abs_type l = (node->getChild(k))->getState();
            temp.push_back(l);
          }
        }
        states_to_explore = temp;
        temp.clear();
      }

      /* explore all the states which are in the outer-approximation
       * of the target but not in the inner-approximation of the target */
      for (size_t p = 0; p < ab; p++) {
        std::vector<abs_type> maybe;
        std::vector<abs_type> maybe_finer;
        for (size_t i = 0; i < Xs_[p]->size(); i++) {
          if (tree_->getMarkingStatus(p, i) == 1) {
            maybe.push_back(i);
          }
        }
        // std::vector<abs_type> temp;
        tree_->finer(p, ab, maybe, maybe_finer);
        for (size_t i = 0; i < maybe_finer.size(); i++) {
          if (tree_->getMarkingStatus(ab, maybe_finer[i])<2) {
            states_to_explore.push_back(maybe_finer[i]);
          }
        }
      }

      // debug purpose: print the intermediate states_to_explore states
      Goal* G = new Goal(*Xs_[ab],states_to_explore);
      intermediate_states_to_explore.push_back(G);
      // end


      Abstraction<X_type, U_type> abs(*Xs_[ab], *U_);
      // cout << "Starting computation of transition relation of layer " << ab << ": \n";
      abs.compute_gb_partial(*Ts_[ab], states_to_explore, sysNext, radNext, *solvers_[ab], avoid);

      // const abs_ptr_type ntr = Ts_[ab]->get_no_transitions();
      // if (ntr!=0) {
      //   double** tr_domain = new double*[ntr];
      //   for(int i = 0; i < ntr; ++i)
      //       tr_domain[i] = new double[2*(*system_->dimX_)+ (*system_->dimU_)];
      //
      //   Ts_[ab]->get_domain(Xs_[ab], U_, tr_domain);
      // }
      // delete[] val;
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
    // BDD uPre(BDD Z, int ab) {
    //     BDD Z2 = Z.Permute(permutesXtoX2_[0]);
    //     BDD uPreZ = uTs_[ab]->symbolicSet_.AndAbstract(Z2, *notXUvars_[0]);
    //     uPreZ = uPreZ.ExistAbstract(*notXvars_[0]);
    //     return uPreZ;
    // }


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
    void initialize(System* system, isG G, isO O, bool verbose=false) {
        //ddmgr_ = new Cudd;
        system_ = system;
        verbose_ = verbose;

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
        verifySave(G);

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
          tree_->markNode(numAbs-1, i);
          // debug purpose
          // std::array<double, 2> x;
          // Xs_[numAbs-1]->itox(i, x);
          // std::cout << "goal: " << i << " : " << x[0] << ", " << x[1] << '\n';
        }
      }
      saveGO();
      // auto target = [&](const abs_type &abs_state, const int ab) {
      //   if (tree_->getMarkingStatus(ab, abs_state))
      //     return true;
      //   else
      //     return false;
      // };
    }

    /*! Initializes SymbolicSet data members.
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     *  \param[in]	addG	Function pointer specifying the points that should be added to the goal set.
     */
    void initializeSymbolicSets() {
      checkMakeDir("X");
        for (int i = 0; i < *system_->numAbs_; i++) {
    			UniformGrid* X = new UniformGrid(*system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[i]);
    			Xs_.push_back(X);
    			//X2s_.push_back(X);
          std::string file = "X/X" + std::to_string(i);
          write_to_file(*X, file);
          // debug purpose
          // std::cout << "Layer " << i << "state space: " << '\n';
          // X->print_info();
        }
		clog << "Xs_ and X2s_ initialized with full domain.\n";

      tree_ = new StateTree(*system_->numAbs_, Xs_, system_->etaRatio_);

        U_ = new UniformGrid(*system_->dimU_, system_->lbU_, system_->ubU_, system_->etaU_);
        write_to_file(*U_, "U");
        clog << "U_ initialized with full domain.\n";

        Xs_[0]->print_info();

		//delete U;

      if (verbose_) {
        for (int i=0; i<*system_->numAbs_; i++) {
          std::cout << "State space of layer " << i << ":" << '\n';
          Xs_[i]->print_info();
        }
        std::cout << "Input space:" << '\n';
        U_->print_info();
      }

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
      if (i!=0) {
        Ts->init_infrastructure(Xs_[i]->size(),U_->size());
        Ts->init_transitions(0);
      }
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

    /*! Saves and prints to log file some information related to the reachability/always-eventually specification. */
    template<class classG>
    void verifySave(classG G) {
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
        write_to_file(*Xs_[0], "plotting/X");
        write_to_file(*Xs_[*system_->numAbs_-1], avoid, "plotting/O");
        write_to_file(*Xs_[*system_->numAbs_-1], G, "plotting/G");
        // checkMakeDir("O");
        // saveVec(Os_, "O/O");
        // checkMakeDir("G");
        // saveVec(Gs_, "G/G");
    }

    /* copied from GameSolver.hh of SCOTSv0.2 (while suitably adapting the functions target and obstacle)
     * @brief solve reachability game according to Algorithm 2 in  <a href="./../../manual/manual.pdf">manual</a>
     *
     * @param[in] trans_function - TransitionFunction of the symbolic model
     * @param[in] target - lambda expression of the form
     *                      \verbatim [] (const abs_type &i) -> bool \endverbatim
     *                      returns true if state i is in target set and false otherwise
     *
     * @param[in] avoid  - OPTIONALLY provide lambda expression of the form
     *                      \verbatim [] (const abs_type &i) -> bool \endverbatim
     *                      returns true if state i is in avoid set and false otherwise
     *
     * @param[out] value - OPTIONALLY provide std::vector<double> value to obtain the value function
     *
     * @return -  WinningDomain that contains the set of winning states and valid inputs
     **/
    // template<class F1, class F2=decltype(params::avoid)>
    WinningDomain solve_reachability_game(int ab,
                                          std::vector<double> & value = params::value ) {
      TransitionFunction& trans_function = *Ts_[ab];
      const UniformGrid* ss = Xs_[ab];
      /* size of state alphabet */
      abs_type N=trans_function.m_no_states;
      /* size of input alphabet */
      abs_type M=trans_function.m_no_inputs;

      /* used to encode that a state is not in the winning domain */
      abs_type loosing = std::numeric_limits<abs_type>::max();
      if(M > loosing-1) {
        throw std::runtime_error("scots::solve_reachability_game: Number of inputs exceeds maximum supported value");
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
        if(target(i,ab) && !avoid(i,ss)) {
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

      /* main loop */
      while(!fifo.empty()) {
        /* get state to be processed */
        abs_type q=fifo.front();
        fifo.pop();
        /* loop over each input */
        for(abs_type j=0; j<M; j++) {
          /* loop over pre's associated with this input */
          for(abs_ptr_type v=0; v<trans_function.m_no_pre[q*M+j]; v++) {
          	std::vector<abs_type> debug = trans_function.get_pre(q, j);
          	abs_type i=trans_function.m_pre[trans_function.m_pre_ptr[q*M+j]+v];
            if(avoid(i,ss))
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

      /* if the default value function was used, free the memory of the static object*/
      if(value == scots::params::value){
          value.clear();
          value.shrink_to_fit();
      }

      return WinningDomain(N,M,std::move(win_domain),std::vector<bool>{},loosing);
    }


    /*! Saves a snapshot of a controller and its domain into the sequence of final controllers and controller domains.
       \param[in] ab	0-index of the abstraction which the controller and controller domain that should be saved belong to.
       \param[in] w   1-winning domain of the controller.
    */
    void saveCZ(const int ab, WinningDomain& w) {
       std::vector<abs_type>* Z = new std::vector<abs_type>();
       // if (!tree_->getMarkedStates(ab, *Z))
       //  return;
       // Goal* Zg = new Goal(*Xs_[ab], *Z);

       // debug purpose
       if (!tree_->getMarkedStates(*system_->numAbs_-1, *Z))
        return;
       Goal* Zg = new Goal(*Xs_[*system_->numAbs_-1], *Z);
       // end

       finalZs_.push_back(Zg);

       StaticController* C = new StaticController(*Xs_[ab], *U_, std::move(w), *tau_[ab]);
       finalCs_.push_back(C);
       finalAbs_.push_back(ab);
    }

    /*! Saves the goal and obstacle states
       \param[in] ab	0-index of the abstraction which the controller and controller domain that should be saved belong to.
       \param[in] w   1-winning domain of the controller.
    */
    void saveGO() {
      for (size_t i = 0; i < *system_->numAbs_; i++) {
        std::vector<abs_type>* G = new std::vector<abs_type>();
        bool flag = tree_->getMarkedStates(i, *G);

        Goal* Gg = new Goal(*Xs_[i], *G);
        Gs_.push_back(Gg);

        //debug purpose
        // for (size_t j = 0; j < Gg->points.size(); j++) {
        //   std::cout << "marked goal: " << Gg->points[j] << '\n';
        // }
        // std::cout << "debug\n\n";
        // Gg->state_grid.print_info();
        // abs_type id = 57;
        // std::vector<double> x;
        // Gg->state_grid.itox(id,x);
      }
       //
       //
       // StaticController* C = new StaticController(*Xs_[ab], *U_, std::move(w));
       // finalCs_.push_back(C);
       // finalAbs_.push_back(ab);
    }
};
}

#endif /* ADAPTABSREACH_HH_ */
