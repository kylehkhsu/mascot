/*! \file System.hh
 *  Contains the System class. */

#ifndef SYSTEM_HH_
#define SYSTEM_HH_

using namespace helper;

namespace scots {

/*! \class System
 *  \brief A class used for storing and passing abstraction parameters to an Adaptive object.
 */
class System {
public:
    int* dimX_; /*!< Dimensionality of state space. */
    double* lbX_; /*!< Lowermost grid point of state space. */
    double* ubX_; /*!< Uppermost grid point of state space. */
    double* etaX_; /*!< Coarsest grid spacing of state space. */
    double* tau_; /*!< Longest sampling time step. */
    int* dimU_; /*!< Dimensionality of input space. */
    double* lbU_; /*!< Lowermost grid point of input space. */
    double* ubU_; /*!< Uppermost grid point of input space. */
    double* etaU_; /*!< Grid spacing of input space abstraction in each dimension. */

    /*! Constructor for a System object.
        \param[in]	dimX		Dimensionality of the state space.
        \param[in]	lbX			Lowermost grid point of the state space.
        \param[in]	ubX 		Uppermost grid point of the state space.
        \param[in]	etaX 		Coarsest grid spacing of the state space.
        \param[in]	tau 		Coarsest time step.
        \param[in]	dimU		Dimensionality of the input space.
        \param[in]	lbU 		Lowermost grid point of the input space.
        \param[in]	ubU			Uppermost grid point of the input space.
        \param[in] 	etaU		Grid spacing of the input space.
    */
    System(int dimX, double* lbX, double* ubX, double* etaX, double tau,
          int dimU, double* lbU, double* ubU, double* etaU) {
        dimX_ = new int;
        lbX_ = new double[dimX];
        ubX_ = new double[dimX];
        etaX_ = new double[dimX];
        tau_ = new double;
        dimU_ = new int;
        lbU_ = new double[dimU];
        ubU_ = new double[dimU];
        etaU_ = new double[dimU];

        *dimX_ = dimX;
        for (int i = 0; i < dimX; i++) {
            lbX_[i] = lbX[i];
            ubX_[i] = ubX[i];
            etaX_[i] = etaX[i];
        }
        *tau_ = tau;
        *dimU_ = dimU;
        for (int i = 0; i < dimU; i++) {
            lbU_[i] = lbU[i];
            ubU_[i] = ubU[i];
            etaU_[i] = etaU[i];
        }
    }

    /*! Destructor for a System object. */
    ~System() {
        delete dimX_;
        delete[] lbX_;
        delete[] ubX_;
        delete[] etaX_;
        delete tau_;
        delete dimU_;
        delete[] lbU_;
        delete[] ubU_;
        delete[] etaU_;
    }
};
}

#endif /* SYSTEM_HH_ */
