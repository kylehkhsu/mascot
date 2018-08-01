/*
 * scots.hh
 *
 *     created: Jan 2017
 *      author: Matthias Rungger
 */

/** @file **/

#ifndef SCOTS_HH_
#define SCOTS_HH_

#include "TransitionFunction.hh"
#include "UniformGrid.hh"
#include "Abstraction.hh"
//#include "GameSolver.hh"
#include "WinningDomain.hh"
//#include "StaticController.hh"
#include "InputOutput.hh"
#include "AdaptAbsReach.hh"
//#include "Helper.hh"
#include "Goal.hh"


/* if scots is used in connection with the cudd library */
#ifdef  SCOTS_BDD
/* cudd library */
#include "dddmp.h"
/* scots classes with bdd support */
#include "SymbolicSet.hh"
#include "SymbolicModel.hh"
#include "EnfPre.hh"
#endif

#endif /* SCOTS_HH_ */

