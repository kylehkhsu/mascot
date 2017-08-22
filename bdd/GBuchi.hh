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






};

}

#endif /* GBUCHI_HH_ */
