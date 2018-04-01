/*! \file System.hh
 *  Contains the System class. */

#ifndef SYSTEM_HH_
#define SYSTEM_HH_

#include "Helper.hh"

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
    double* etaRatio_; /*!< Ratio between state space grid spacings of consecutive abstractions. */
    double* tauRatio_; /*!< Ratio between time steps of consecutive abstractions. */
    int* nSubInt_; /*!< Number of sub-intervals in ODE solving per time step. */
    int* numAbs_; /*!< Number of abstractions of different granularity. */

    /*! Constructor for a System object.
     *  \param[in]	dimX		Dimensionality of the state space.
     *  \param[in]	lbX			Lowermost grid point of the state space.
     *  \param[in]	ubX 		Uppermost grid point of the state space.
     *  \param[in]	etaX 		Coarsest grid spacing of the state space.
     *  \param[in]	tau 		Coarsest time step.
     *  \param[in]	dimU		Dimensionality of the input space.
     *  \param[in]	lbU 		Lowermost grid point of the input space.
     *  \param[in]	ubU			Uppermost grid point of the input space.
     *  \param[in] 	etaU		Grid spacing of the input space.
     *  \param[in]	etaRatio	Ratio between grid spacings of the input space in consecutive abstractions.
     *  \param[in]	tauRatio	Ratio between time steps of consecutive abstractions.
     *  \param[in]  nSubInt     Number of sub-intervals in ODE solving per time step.
     *  \param[in]	numAbs		Number of abstractions.
     */
    System(int dimX, double* lbX, double* ubX, double* etaX, double tau,
           int dimU, double* lbU, double* ubU, double* etaU,
           double* etaRatio, double tauRatio, int nSubInt, int numAbs) {
        dimX_ = new int;
        lbX_ = new double[dimX];
        ubX_ = new double[dimX];
        etaX_ = new double[dimX];
        tau_ = new double;
        dimU_ = new int;
        lbU_ = new double[dimU];
        ubU_ = new double[dimU];
        etaU_ = new double[dimU];
        etaRatio_ = new double[dimX];
        tauRatio_ = new double;
        nSubInt_ = new int;
        numAbs_ = new int;

        *dimX_ = dimX;
        for (int i = 0; i < dimX; i++) {
            lbX_[i] = lbX[i];
            ubX_[i] = ubX[i];
            etaX_[i] = etaX[i];
            etaRatio_[i] = etaRatio[i];
        }
        *tau_ = tau;
        *dimU_ = dimU;
        for (int i = 0; i < dimU; i++) {
            lbU_[i] = lbU[i];
            ubU_[i] = ubU[i];
            etaU_[i] = etaU[i];
        }
        *tauRatio_ = tauRatio;
        *nSubInt_ = nSubInt;
        *numAbs_ = numAbs;
    }
    System(int dimX, double* lbX, double* ubX, double* etaX, double tau,
           double* etaRatio) {
        dimX_ = new int;
        lbX_ = new double[dimX];
        ubX_ = new double[dimX];
        etaX_ = new double[dimX];
        tau_ = new double;
        dimU_ = new int;
        *dimU_ = 1;
        lbU_ = new double[*dimU_];
        ubU_ = new double[*dimU_];
        etaU_ = new double[*dimU_];
        etaRatio_ = new double[dimX];
        tauRatio_ = new double;
        nSubInt_ = new int;
        numAbs_ = new int;

        *dimX_ = dimX;
        for (int i = 0; i < dimX; i++) {
            lbX_[i] = lbX[i];
            ubX_[i] = ubX[i];
            etaX_[i] = etaX[i];
            etaRatio_[i] = etaRatio[i];
        }
        *tau_ = tau;
        *lbU_ = 0;
        *ubU_ = 1;
        *etaU_ = 1;
        *tauRatio_ = -1;
        *nSubInt_ = -1;
        *numAbs_ = -1;
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
        delete[] etaRatio_;
        delete tauRatio_;
        delete nSubInt_;
        delete numAbs_;
    }
};
}

#endif /* SYSTEM_HH_ */
