/**
 * @file
 * This file contains the implementation of BinaryState, which is
 * the base class for the Standard character data type in RevBayes.
 *
 * @brief Implementation of BinaryState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-05-24 09:58:04 +0200 (Thu, 24 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: BinaryState.cpp 1568 2012-05-24 07:58:04Z hoehna $
 */

#include "BinaryState.h"

using namespace RevBayesCore;

/** Default constructor */
BinaryState::BinaryState(size_t n) : StandardState("0", "01")
{

}


/** Default constructor */
BinaryState::BinaryState(const std::string &s) : StandardState(s, "01")
{ 
}


