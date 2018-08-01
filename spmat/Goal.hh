/*
 * Goal.hh
 *
 * created: July, 2018
 *  author: Kaushik Mallik
 */

 /** @file **/
 #ifndef GOAL_HH_
 #define GOAL_HH_

 // #include <iostream>
 // #include <cstring>
 // #include <memory>
 #include <vector>
 #include <string> // debug purpose

 /* to get abs_type alias */
 #include "Goal.hh"
 #include "UniformGrid.hh"

 // #include "mex.h" // debug purpose

 // #include "StateTreeNode.hh"

 // using std::clog;
 // using std::freopen;
 // using std::string;
 // using std::vector;

  /** @namespace scots **/
  namespace scots {
    /**
     * @class Goal
     *
     * @brief Intermediate goals during multi-layered synthesis
     *
     **/
     class Goal {
     public:
       /* full state space */
       UniformGrid state_grid;
       /* points */
       std::vector<abs_type> points;
     public:
       /* default constructor */
       Goal() : state_grid(UniformGrid()),
                points({}) {}
        /* destructor */
        virtual
        ~Goal() = default;
        /* copy constructor */
        Goal(const Goal& other) : Goal() {
          *this=other;
        }
        /* move constructor */
        Goal(Goal&& other) : Goal() {
          *this=std::move(other);
        }
        /* copy assignment operator */
        Goal& operator=(const Goal &other) {
          if(this==&other)
            return *this;
          state_grid = UniformGrid(other.state_grid);
          if(state_grid.get_dim() != 0) {
            points.clear();
            for(int i=0; i<state_grid.get_dim(); i++) {
              points.push_back(other.points[i]);
            }
          }
          return *this;
        }
        /* move assignment operator */
        Goal& operator=(Goal&& other) {
          state_grid=std::move(other.state_grid);
          points=std::move(other.points);
          return *this;
        }

       /* constructor */
       Goal(const UniformGrid& grid,
            std::vector<abs_type> array) :
            state_grid(grid),
            points(array) {state_grid.print_info();}

       // std::vector<double> get_eta() const {
       //   std::vector<double> eta;
       //   for(int i=0; i<m_dim; i++) {
       //     eta.push_back(m_eta[i]);
       //   }
       //   return eta;
       // }

       /* get a std::vector containing the states in the target **/
       template<class state_type>
       std::vector<state_type> get_domain() const {
         /* convert points to cell centers */
         std::vector<std::vector<double>> domain(points.size());
         for(abs_type i=0; i<points.size(); i++) {
           state_type x;
           state_grid.itox(points[i],x);
           domain[i]=x;
           // debug purpose
           // std::string Str = "domain: ";
           // Str += std::to_string(points[i]);
           // Str += ": ";
           // Str += std::to_string(x[0]);
           // Str += ", ";
           // Str += std::to_string(x[1]);
           // Str += "\n";
           // int n = Str.length();
           // char char_array[n+1];
           // strcpy(char_array, Str.c_str());
           // int h = mexPrintf(char_array);
           // state_grid.print_info();
           // end
         }
         return domain;
       }

      /* check if a given (continuous) state is in the goal or not */
      template<class point>
      bool check_state(const point& x) const {
        abs_type i = state_grid.xtoi(x);
        for (size_t j = 0; j < points.size(); j++) {
          // if (i==state_grid.xtoi(points[j])) {
          if (i==points[j]) {
            return true;
          }
        }
        return false;
      }
     };
   }

   #endif
