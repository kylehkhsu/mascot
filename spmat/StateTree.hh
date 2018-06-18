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
       // std::array<std::array<StateTreeNode*>*> db_;
       vector<StateTreeNode**> db_;
     public:
       /* constructor */
       StateTree(const int numAbs,
                 const std::vector<UniformGrid*> ss,
                 const double* eta_ratio) :
              numAbs_(numAbs) {
                /* number of cells in different layers */
                abs_type no_states_[numAbs_];
                for (int ab = 0; ab < numAbs_; ab++)
                  no_states_[ab] = ss[ab]->size();

                  /* debug purpose */
                  std::cout << "no. of states in different layers: " << no_states_[0] << ", " << no_states_[1] << ", " << no_states_[2] << '\n';
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
                  // db_[ab] = nodePtrArray;
                  db_.push_back(nodePtrArray);
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
                for (int ab = 0; ab < numAbs_-1; ab++) {
                  // std::vector<abs_type> no_grid_points_ab = ss[ab].get_no_gp_per_dim();
                  // std::vector<abs_type> no_grid_points_nextAb = ss[ab+1].get_no_gp_per_dim();
                  abs_type head = 0; // next state in the next layer which is to be assigned
                  for (abs_type i = 0; i < no_states_[ab]; i++) {
                    // if (ab == 0) {
                    //   db_[ab][i]->parent_ = nullptr;
                    // }
                    // elseif (ab == numAbs_-1) {
                    //   db_[ab][i]->no_child = 0;
                    // }
                    // else {
                    db_[ab][i]->no_child_ = no_child;
                    for (int j = 0; j < no_child; j++) {
                      db_[ab][i]->child_.push_back(db_[ab+1][head+j]);
                      db_[ab+1][head+j]->parent_ = db_[ab][i];
                    }
                    head += no_child;
                    // }
                  }
                }
              }
            // }

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
            updateParent(curAb, node, 2);
            updateChildren(curAb, node);
          }

          /* recursively update the ancestors */
          void  updateParent(int curAb, StateTreeNode* node, int val) {
            if (curAb != 0) {
              if (val==2) {
                (node->parent_)->no_m_child_ += 1;
                if ((node->parent_)->no_m_child_ = (node->parent_)->no_child_) {
                  (node->parent_)->marking_ = 2;
                }
                else {
                  (node->parent_)->marking_ = 1;
                }
              }
              else {
                (node->parent_)->marking_ = 1;
              }
              updateParent(curAb-1, node->parent_, (node->parent_)->marking_);
            }
            else return;
          }

          /* recursively update the descendants */
          void  updateChildren(int curAb, StateTreeNode* node) {
            if (curAb != numAbs_) {
              for (int n = 0; n < node->no_child_; n++) {
                (node->child_[n])->marking_ = 2;
                updateChildren(curAb+1, node->child_[n]);
              }
            }
            else return;
          }

          /* return the corresponding node in the tree */
          StateTreeNode* getNode(const int ab, const abs_type state) {
            return db_[ab][state];
          }

          /* return the marking status of a state in the tree */
          int getMarkingStatus(const int ab, const abs_type state) {
            return db_[ab][state]->marking_;
          }

     }; /* close class definitioin */
  } /* close namespace */

  #endif /* STATETREE_HH_ */
