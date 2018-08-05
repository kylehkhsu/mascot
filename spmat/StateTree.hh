/*
 * StateTree.hh
 *
 * created: May 4, 2018
 *  author: Kaushik Mallik
 */

 /** @file **/
 #ifndef STATETREE_HH_
 #define STATETREE_HH_

 #include <iostream>
 #include <cstring>
 #include <memory>
 #include <vector>
 #include "TicToc.hh"

 /* to get abs_type alias */
 #include "UniformGrid.hh"

 #include "StateTreeNode.hh"

 using std::clog;
 using std::freopen;
 using std::string;
 using std::vector;

  /** @namespace scots **/
  namespace scots {
    /**
     * @class StateTree
     *
     * @brief Arrangement of StateTreeNodes of different layers in the form of a tree for fast mapping
     *
     **/
     class StateTree {
     protected:
       /* number of layers */
       int numAbs_;
       /* array containing number of abstract states/cells in different layers */
       // std::array<abs_type> no_states_;
       abs_type* no_states_;
       /* array[numAbs_] containing pointers to arrays containing pointers to StateTreeNodes */
       StateTreeNode*** db_;
       // vector<StateTreeNode**> db_;
     public:
       /* constructor */
       StateTree(const int numAbs,
                 const std::vector<UniformGrid*> ss,
                 const double* eta_ratio) :
              numAbs_(numAbs) {
                /* number of cells in different layers */
                no_states_ = new abs_type[numAbs_];
                for (int ab = 0; ab < numAbs_; ab++)
                  no_states_[ab] = ss[ab]->size();

                /* initialize db_ */
                db_ = new StateTreeNode**[numAbs_];
                for (size_t i = 0; i < numAbs_; i++) {
                  db_[i] = new StateTreeNode*[no_states_[i]];
                }
                /* debug purpose */
                // for (int i=0; i<numAbs_; i++)
                //   std::cout << "no. of states in layer " << i  << ": "<< no_states_[i] << '\n';
                /* end */
                /* first pass: initialization of tree nodes */
                StateTreeNode** nodePtrArray = nullptr;
                for (int ab = 0; ab < numAbs_; ab++) {
                  // StateTreeNode* nodePtrArray[no_states_[ab]];
                  nodePtrArray = new StateTreeNode*[no_states_[ab]];
                  for (abs_type i = 0; i < no_states_[ab]; i++) {
                    // StateTreeNode node(ab, i);
                    // nodePtrArray[i] = &node;
                    StateTreeNode* node = new StateTreeNode(ab, i);
                    nodePtrArray[i] = node;
                  }
                db_[ab] = nodePtrArray;
                // db_.push_back(nodePtrArray);
                }
                /* second pass: mapping of tree nodes across layers */
                int dim = ss[0]->get_dim(); /* the dimension is same in all layers */
                int no_child = 1;
                for (int d = 0; d < dim; d++) {
                  // for (int j = 0; j < eta_ratio[d]; j++)
                  no_child = no_child * eta_ratio[d];
                }
                /* loop over layers 0 to numAbs_ -1
                 * map each state of present state to no_child of consecutive
                 * states of the next layer */
                 abs_type quo;
                 abs_type first_child;
                 abs_type child_id;
                 for (int ab = 0; ab < numAbs_-1; ab++) {
                   std::queue<std::vector<int>> fifo;
                   for (size_t i = 0; i < dim; i++) {
                     std::vector<int> v;
                     int p = 0;
                     for (size_t j = 0; j < eta_ratio[i]; j++) {
                       v.push_back(p);
                       p += ss[ab+1]->m_NN[i];
                     }
                     fifo.push(v);
                   }
                   fifo = get_offset(fifo);
                   std::vector<int> offset = fifo.front();
                   // std::cout << "Abstraction: " << ab << "\n\n";
                   for (abs_type i = 0; i < no_states_[ab]; i++) {
                     // std::cout << "State: " << i << " : " << '\n';
                     db_[ab][i]->no_child_ = no_child;
                     abs_type id = i;
                     first_child = 0;
                     for (int j = dim-1; j > 0; j--) {
                       quo = id/ss[ab]->m_NN[j];
                       first_child += quo*eta_ratio[j]*ss[ab+1]->m_NN[j];
                       id = id%ss[ab]->m_NN[j];
                     }
                     first_child += id*eta_ratio[0];
                     // std::cout << "Children: " << '\n';
                     for (size_t j = 0; j < offset.size(); j++) {
                       child_id = first_child + offset[j];
                       // std::cout << "\t " << child_id << '\n';
                       db_[ab][i]->child_.push_back(db_[ab+1][child_id]);
                       db_[ab+1][child_id]->parent_ = db_[ab][i];
                     }
                     // for (int j = 0; j < dim; j++) {
                     //   for (int k = 0; k < eta_ratio[j]; k++) {
                     //     child_id = first_child + k;
                     //     std::cout << "\t " << child_id << '\n';
                     //     db_[ab][i]->child_.push_back(db_[ab+1][child_id]);
                     //     db_[ab+1][child_id]->parent_ = db_[ab][i];
                     //   }
                     //   if (j < dim-1) {
                     //    first_child += ss[ab+1]->m_NN[j+1];
                     //   }
                     // }

                   }
                 }
                // for (int ab = 0; ab < numAbs_-1; ab++) {
                //   abs_type head = 0; // next state in the next layer which is to be assigned
                //   for (abs_type i = 0; i < no_states_[ab]; i++) {
                //     db_[ab][i]->no_child_ = no_child;
                //     for (int j = 0; j < no_child; j++) {
                //         db_[ab][i]->child_.push_back(db_[ab+1][head+j]);
                //         db_[ab+1][head+j]->parent_ = db_[ab][i];
                //     }
                //     head += no_child;
                //   }
                // }
              }

    // /* copy constructor */
    // StateTree(const StateTree& other) : StateTree() {
    //   *this = other;
    // }

    /* copy assignment operator */
    StateTree& operator=(const StateTree &other) {
      numAbs_ = other.numAbs_;
      no_states_ = other.no_states_;
      db_ = other.db_;

      return *this;
    }

          /** @brief update the marking of the tree **/
          void markNode(const int curAb, const abs_type state) {
            /* inputs: the node <node> in layer <curAb> which was marked with 2
             * by the fixpoint algorithm */
            /* action: mark all the descendants with 2
             * mark all the ancestors with either 2 or 1, depending on whether
             * all the children of one particular ancestor are marked 2 or not, respectively */
            StateTreeNode* node = getNode(curAb, state);
            node->marking_ = 2;
            updateParent(node, node->marking_);
            updateChildren(node);
          }

          /* recursively update the ancestors
           * node is the child with marking = val, whose parent has to be updated
           * the input parameter val can never be 0 */
          // void  updateParent(int curAb, StateTreeNode* node, int val) {
          void  updateParent(StateTreeNode* node, int val) {
            // if (curAb != 0) {
            if (node->parent_ != nullptr) {
              int oldMarking = (node->parent_)->marking_;
              if (val==2) {
                (node->parent_)->no_m_child_ += 1;
                if ((node->parent_)->no_m_child_ == (node->parent_)->no_child_) {
                  (node->parent_)->marking_ = 2;
                }
                else {
                  (node->parent_)->marking_ = 1;
                }
              }
              else { // val = 1
                (node->parent_)->marking_ = 1;
              }
              if (oldMarking != (node->parent_)->marking_) { // parent changed marking status
                updateParent(node->parent_, (node->parent_)->marking_);
              }
            }
            return;
          }

          /* recursively update the descendants */
          void  updateChildren(StateTreeNode* node) {
            node->no_m_child_ = node->no_child_;
            if (node->no_child_ != 0) {
              for (int n = 0; n < node->no_child_; n++) {
                int oldMarking = (node->child_[n])->marking_;
                (node->child_[n])->marking_ = 2;
                if (oldMarking != (node->child_[n])->marking_)
                  updateChildren(node->child_[n]);
              }
            }
            return;
          }

          /* return the corresponding node in the tree */
          StateTreeNode* getNode(const int ab, const abs_type state) {
            return db_[ab][state];
          }

          /* return the marking status of a state in the tree */
          int getMarkingStatus(const int ab, const abs_type state) {
            if (state >= no_states_[ab]) {
              std::ostringstream os;
              os << "\nscots::StateTree: the state " << state << " is outside the state space of layer " << ab ;
              throw std::runtime_error(os.str().c_str());
            }
            return db_[ab][state]->marking_;
          }

          /* construct a vector of states in a given layer which are marked with 2 */
          bool getMarkedStates(const int ab, std::vector<abs_type>& marked) {
            bool flag = false; // flag = false iff no state is marked with 2
            for (size_t i = 0; i < no_states_[ab]; i++) {
              if (getMarkingStatus(ab, i)==2) {
                marked.push_back(i);
                flag = true;
                // debug purpose
                // std::cout << "marked state: " << i << '\n';
              }
            }
            return flag;
          }
          /* downward projection between two consecutive layers */
          void nextFiner(const int curAb, std::vector<abs_type>& sc, std::vector<abs_type>& sf) {
            StateTreeNode *q, *ch;
            for (size_t i = 0; i < sc.size(); i++) {
              q = db_[curAb][sc[i]];
              for (int j=0; j<q->get_no_child(); j++) {
                ch = q->child_[j];
                sf.push_back(ch->getState());
              }
            }
          }
          /* downward projection between any two layers */
          void finer(const int curAb, const int targetAb, std::vector<abs_type>& sc, std::vector<abs_type>& sf) {
            for (size_t i = curAb; i < targetAb; i++) {
              sf.clear();
              nextFiner(i, sc, sf);
              sc = sf;
            }
          }
          /* prints some information */
          void print_info() {
            std::cout << "\n State Tree:" << '\n';
            std::cout << "No. of layers:  " << numAbs_ << '\n';
            std::cout << "No. of states in each layer:  ";
            for (size_t i = 0; i < numAbs_; i++) {
              std::cout << no_states_[i] << " ";
            }
            std::cout << '\n';
          }

        private:
          /* compute offset vector */
          template<class T>
          std::queue<std::vector<T>> get_offset(std::queue<std::vector<T>> fifo) {
            if (fifo.size()>2) {
              std::queue<std::vector<T>> temp;
              temp.push(fifo.front());
              fifo.pop();
              temp.push(fifo.front());
              fifo.pop();
              std::queue<std::vector<T>> merged = get_offset(temp);
              fifo.push(merged.front());
              fifo = get_offset(fifo);
            }
            if (fifo.size()==2) {
              std::vector<T> v1 = fifo.front();
              fifo.pop();
              std::vector<T> v2 = fifo.front();
              fifo.pop();
              std::vector<T> vmerged;
              for (size_t i = 0; i < v1.size(); i++) {
                for (size_t j = 0; j < v2.size(); j++) {
                  vmerged.push_back(v1[i]+v2[j]);
                }
              }
              fifo.push(vmerged);
              return fifo;
            }
            return fifo;
          }

     }; /* close class definitioin */
  } /* close namespace */

  #endif /* STATETREE_HH_ */
