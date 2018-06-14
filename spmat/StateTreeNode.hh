/*
 * StateTreeNode.hh
 *
 * created: May 4, 2018
 *  author: Kaushik Mallik
 */

 /** @file **/
 #ifndef STATETREENODE_HH_
 #define STATETREENODE_HH_

 #include <iostream>
 #include <cstring>
 #include <memory>
 #include <vector>
 #include "TicToc.hh"

 using std::clog;
 using std::freopen;
 using std::string;
 using std::vector;

 /* to get abs_type alias */
 #include "UniformGrid.hh"

 /** @namespace scots **/
 namespace scots {
   /**
    * @class StateTreeNode
    *
    * @brief Encapsulation of inheritence information of states
    *
    **/
    class StateTreeNode {
    private:
      /* layer */
      const int layer_;
      /* state */
      const abs_type state_;
      /* number of children, 0 if leaf */
      int no_child_;
      /* pointer to children, null if leaf */
      vector<StateTreeNode*> child_;
      /* pointer to parent, null if root */
      StateTreeNode* parent_;
      /* marking status
       *  2: fully marked, which causes all children to be fully marked (recursively),
       *  1: some children have marking >= 1, and there is at least one child with marking <= 1,
       *  0: otherwise */
      int marking_;
      /* number of children marked 2; for non-leaf nodes, whenever no_m_child = no_child, marking_ goes to 2 */
      int no_m_child_;
      friend class StateTree;
    public:
      /* @cond EXCLUDE from doxygen */
      /* default constructor (must contain a state) */
      StateTreeNode(const int layer, const abs_type state) :  layer_(layer),
                                                      state_(state),
                                                      no_child_(0),
                                                      parent_(nullptr),
                                                      marking_(0),
                                                      no_m_child_(0) {}
      // /* move constructor */
      // StateTreeNode(StateTreeNode&& other) {
      //   *this = std::move(other);
      // }
      // /* move assignment operator */
      // StateTreeNode& operator = (StateTreeNode&& other) {
      //   layer_ = std::move(other.layer_);
      //   state_ = std::move(other.state_);
      //   no_child_ = std::move(other.no_child_);
      //   child_ = std::move(other.child_);
      //   parent_ = std::move(other.parent_);
      //   marking_ = std::move(other.marking_);
      //   no_m_child_ = std::move(other.no_m_child_);
      //
      //   return *this;
      // }
      /* Destructor */
      ~StateTreeNode() {
        // delete layer_;
        // delete state_;
        // delete no_child;
        deleteVec(child_);
        delete parent_;
        // delete marking_;
        // delete no_m_child;
      }
      /* @endcond */
      /* constructor */
      StateTreeNode(const int layer,
                    const abs_type state,
                    int no_child,
                    vector<StateTreeNode*> child,
                    StateTreeNode* parent):
              layer_(layer), state_(state), no_child_(no_child), child_(child), parent_(parent), marking_(0), no_m_child_(0) {}
  }; /* end of class definition */
} /* end of namespace */

#endif /* STATETREENODE_HH_ */
