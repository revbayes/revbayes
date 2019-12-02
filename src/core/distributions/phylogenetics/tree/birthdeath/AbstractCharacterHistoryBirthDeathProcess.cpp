//
//  AbstractCharacterHistoryBirthDeathProcess.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/10/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "AbstractCharacterHistoryBirthDeathProcess.h"
//#include "Clade.h"
//#include "RandomNumberFactory.h"
//#include "RandomNumberGenerator.h"
//#include "RbConstants.h"
//#include "RbMathCombinatorialFunctions.h"
//#include "StochasticNode.h"
//#include "TopologyNode.h"
//
//#include <algorithm>
//#include <cmath>

using namespace RevBayesCore;

AbstractCharacterHistoryBirthDeathProcess::AbstractCharacterHistoryBirthDeathProcess( void ) : TypedDistribution<Tree>( new Tree() )
{
    ; // do nothing
}


AbstractCharacterHistoryBirthDeathProcess::~AbstractCharacterHistoryBirthDeathProcess()
{
    ; // do nothing
}
