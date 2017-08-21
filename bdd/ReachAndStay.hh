/*! \file ReachAndStay.hh
 *  Contains the ReachAndStay class. */

#ifndef REACHANDSTAY_HH_
#define REACHANDSTAY_HH_

#include "Reach.hh"
#include "Safe.hh"

using namespace helper;

namespace scots {

/*! \class ReachAndStay
 *  \brief A class (derived from base Adaptive) that does adaptive multiscale abstraction-based synthesis for a reach-and-stay specification.
 */
template<class X_type, class U_type>
class ReachAndStay: public Reach<X_type, U_type>, public Safe<X_type, U_type> {
public:

    /*! Constructor for a ReachAndStay object. */
    ReachAndStay(char* logFile)
                 : Reach<X_type, U_type>(logFile),
                   Safe<X_type, U_type>(logFile),
                   Adaptive<X_type, U_type>(logFile)
    {
    }

    /*! Deconstructor for a ReachAndStay object. */
    ~ReachAndStay() {

    }

    /*! Writes a series of controllers that satisfies a reach-and-stay specification.
     *  \param[in]	addG	Function pointer specifying the points that should be added to the goal set.
     *  \param[in]	addI	Function pointer specifying the points that should be added to the initial set.
     *  \param[in]  startAbs            0-index of the abstraction to start with.
     *  \param[in]	minToGoCoarser		Minimum number of growing fixed point iterations needed before an attempt to go to a coarser abstraction.
     *  \param[in]	minToBeValid		Minimum number of growing fixed point iterations needed before a controller is declared valid.
     *  \param[in]  earlyBreak          If 1, the synthesis ends as soon as I meets the domain of C.
     */
    template<class G_type, class I_type>
    void reachAndStay(G_type addG, I_type addI, int startAbs, int minToGoCoarser, int minToBeValid, int earlyBreak) {
        this->initializeSafe(addG);
        TicToc ttWhole;
        TicToc ttSafe;
        ttWhole.tic();

        // removing obstacles from transition relation
        this->removeFromTs(&(this->Os_));

        // solving safety
        ttSafe.tic();
        for (int curAbs = 0; curAbs < *this->system_->numAbs_; curAbs++) {
            this->nu(curAbs);
        }
        clog << "----------------------------------------safe: ";
        ttSafe.toc();

        checkMakeDir("Z");
        saveVec(this->Zs_, "Z/safeZ");
        checkMakeDir("C");
        saveVec(this->Cs_, "C/safeC");

        this->initializeReach(addI, addI);

        cout << "\nDomain of all safety controllers:\n";
        this->Xs_[*this->system_->numAbs_-1]->printInfo(1);

        // Safe initialized all Zs_ to BDD-1, initializeReach was BS for Gs_.
        for (int i = 0; i < *this->system_->numAbs_; i++) {
            this->Gs_[i]->symbolicSet_ = this->ddmgr_->bddZero();
            this->Zs_[i]->symbolicSet_ = this->ddmgr_->bddZero();
            this->Cs_[i]->symbolicSet_ = this->ddmgr_->bddZero();
        }
        clog << "Gs_, Zs_, Cs_ reset to BDD-0.\n";

        // Xs_[numAbs-1] contains combined domains of all safety controllers in the finest abstraction.
        this->Gs_[*this->system_->numAbs_-1]->symbolicSet_ = this->Xs_[*this->system_->numAbs_-1]->symbolicSet_;

        // Xs_ were changed during safety.
        for (int i = 0; i < *this->system_->numAbs_; i++) {
            this->Xs_[i]->addGridPoints();
        }
        clog << "Xs_ reset to full domain.\n";

        this->Gs_[*this->system_->numAbs_-1]->writeToFile("plotting/G.bdd");

        for (int i = *this->system_->numAbs_-1; i > 0; i--) {
            this->innerCoarserAligned(this->Gs_[i-1], this->Gs_[i], i-1);
            this->Gs_[i-1]->printInfo(1);
        }

        clog << "Gs_ set to projections of domain of safety controllers.\n";

        saveVec(this->Gs_, "G/G");

        // Adaptive safety changes Ts_ and TTs_. Reset.
        this->loadTs();

        this->reach(startAbs, minToGoCoarser, minToBeValid, earlyBreak);

        clog << "----------------------------------------reachAndStay: ";
        ttWhole.toc();

    }


};

}

#endif /* REACHANDSTAY_HH_ */
