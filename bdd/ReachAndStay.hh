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


};

}

#endif /* REACHANDSTAY_HH_ */
