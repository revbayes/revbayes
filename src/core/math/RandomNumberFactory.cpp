/**
 * @file
 * This file contains the implementation of RandomNumberFactory, which is
 * used to manage random number generating objects. The class has a pool
 * of random number objects that it can hand off as needed. This singleton
 * class has two seeds it manages: one is a global seed and the other is
 * a so called local seed.
 *
 * @brief Declaration of RandomNumberFactory
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id$
 */

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbOptions.h"
#include <ctime>
#include <climits>
#include <climits>

using namespace RevBayesCore;

/** Default constructor */
RandomNumberFactory::RandomNumberFactory(void)
{

    seedGenerator = new RandomNumberGenerator();
}


/** Destructor */
RandomNumberFactory::~RandomNumberFactory(void) {

    delete seedGenerator;
}


/** Delete a random number object (remove it from the pool too) */
void RandomNumberFactory::deleteRandomNumberGenerator(RandomNumberGenerator* r) {

    allocatedRandomNumbers.erase( r );
    
    delete r;
}
