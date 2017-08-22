/*! \file GBuchi.hh
 *  Contains the GBuchi class. */

#ifndef GBUCHI_HH_
#define GBUCHI_HH_

#include "Composition.hh"

using std::vector;
using std::cout;

namespace scots {

/*! \class GBuchi
 *  \brief A class (derived from base Composition) that does adaptive multiscale abstraction-based synthesis for a generalized Buchi specification.
 */
class GBuchi: public Composition {
public:


    /*! Constructor for a GBuchi object. */
    GBuchi(char* logFile)
        : Composition(logFile)
    {
    }

//    template<class I_type>
//    void auxAddInitial(I_type addI) {

//    }

    ~GBuchi() {

    }



    BDD Cpre (BDD Z, BDD T, BDD TT, int* permuteXtoX2, BDD notXUVar) {
        // swap to post-state variables
        BDD Z2 = Z.Permute(permuteXtoX2);
        // non-winning post-states
        BDD nZ2 = !Z2;
        // {(x,u)} with non-winning post-state
        BDD nN = T.AndAbstract(nZ2, notXUVar);
        // {(x,u)} with no non-winning post-states
        BDD N = !nN;
        // get rid of non-domain information
        N = TT.AndAbstract(N, notXUVar);
        return N;
    }

    void innerFiner(SymbolicSet* Zc, SymbolicSet* Zf, SymbolicSet* XXc, BDD* notXVarf, SymbolicSet* Xf) {
        BDD Q = XXc->symbolicSet_ & Zc->symbolicSet_;
        Zf->symbolicSet_ = Q.ExistAbstract(*notXVarf) & Xf->symbolicSet_;
    }

    int innerCoarser(SymbolicSet* Zc, SymbolicSet* Zf, SymbolicSet* XXc, BDD* notXVarc, SymbolicSet* Xc, int numFiner) {
        SymbolicSet Qcf(*XXc);
        Qcf.symbolicSet_ = XXc->symbolicSet_ & Zf->symbolicSet_;

        SymbolicSet Qc(*Zc);
        Qc.symbolicSet_ = Qcf.symbolicSet_.ExistAbstract(*notXVarc); // & S1
        Qc.symbolicSet_ &= !(Zc->symbolicSet_); /* don't check states that are already in Zc */

        int* QcMintermWhole;
        SymbolicSet Ccf(Qcf);
        BDD result = ddmgr_->bddZero();
        int* QcMinterm = new int[Qc.nvars_];

        for (Qc.begin(); !Qc.done(); Qc.next()) { // iterate over all coarse cells with any corresponding finer cells in Zf
            QcMintermWhole = (int*) Qc.currentMinterm();
            std::copy(QcMintermWhole + Qc.idBddVars_[0], QcMintermWhole + Qc.idBddVars_[0] + Qc.nvars_, QcMinterm);

            BDD coarseCell = Qc.mintermToBDD(QcMinterm) & Xc->symbolicSet_; // a particular coarse cell

            Ccf.symbolicSet_ = Qcf.symbolicSet_ & coarseCell; // corresponding finer cells to the coarse cell

            if ((Ccf.symbolicSet_.CountMinterm(Ccf.nvars_)) == (numFiner)) { // if there's a full set of finer cells
                result |= coarseCell;
            }
        }

        delete[] QcMinterm;

        if (result == ddmgr_->bddZero()) {
            return 0;
        }

        Zc->symbolicSet_ |= result;
        return 1;
    }







};

}

#endif /* GBUCHI_HH_ */
