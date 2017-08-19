/*! \file Product.hh
 *  Contains the Product class. */

#ifndef PRODUCT_HH_
#define PRODUCT_HH_

#include <vector>
#include "System.hh"
#include "SymbolicModelGrowthBound.hh"
#include "RungeKutta4.hh"

using namespace helper;
using std::vector;

namespace scots {

/*! \class Product
 *  \brief A class used for storing and passing abstraction parameters to an Adaptive object.
 */
class Product {
public:
    Cudd* ddmgr_;
    vector<SymbolicSet*> Ts_;
    int components_;
    vector<SymbolicSet*> Tprods_;

    Product(int components) {
        ddmgr_ = new Cudd;
        components_ = components;
    }

    ~Product() {
        deleteVec(Ts_);
        deleteVec(Tprods_);
        delete ddmgr_;
    }

//    void product(System *dyn, System *pred) {
//        BDD T = Ts_[0]->symbolicSet_;
//        SymbolicSet* Td = new SymbolicSet(*Ts_[0]);
//        Tprods_.push_back(Td);

//        for (int i = 1; i < components_; i++) {
//            T &= Ts_[i]->symbolicSet_;
//            SymbolicSet* Tprev = Tprods_.back();
//            SymbolicSet* Tprod = new SymbolicSet(*Tprev, *Ts_[i]);
//            Tprod->symbolicSet_ = T;
//            Tprods_.push_back(Tprod);
//        }
//        Tprods_[components_-1]->printInfo(1);

//        int permut[44] = {0, 1, 2, 3, 4,
//                          5, 6, 7, 8, 9,
//                          38, 39, 40,
//                          41, 42,
//                          19, 20, 21, 22, 23,
//                          24, 25, 26, 27, 28,
//                          10,
//                          11, 12, 13, 14,
//                          15, 16, 17, 18,
//                          43,
//                          29,
//                          30, 31, 32, 33,
//                          34, 35, 36, 37};

//        Cudd mgr;
//        SymbolicSet X(mgr, *dyn->dimX_, dyn->lbX_, dyn->ubX_, dyn->etaX_, )


//    }

//    void permute() {



//    }

    template<class sys_type, class rad_type, class X_type, class U_type>
    void computeAbstraction(System* system, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        SymbolicSet X(*ddmgr_, *system->dimX_, system->lbX_, system->ubX_, system->etaX_, *system->tau_);
        X.addGridPoints();
        SymbolicSet U(*ddmgr_, *system->dimU_, system->lbU_, system->ubU_, system->etaU_, 0);
        U.addGridPoints();
        SymbolicSet X2(X, 1);

        OdeSolver solver(*system->dimX_, 5, *system->tau_);

        SymbolicModelGrowthBound<X_type, U_type> abs (&X, &U, &X2);
        abs.computeTransitionRelation(sysNext, radNext, solver);

        SymbolicSet C(X, U);
        SymbolicSet* T = new SymbolicSet(C, X2);
        T->symbolicSet_ = abs.transitionRelation_;
        Ts_.push_back(T);

        T->printInfo(1);
    }

};




}

#endif /* Product_HH_ */
