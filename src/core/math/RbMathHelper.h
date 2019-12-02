/**
 * @file
 * This file contains commonly used math functions that are used
 * in RevBayes. 
 *
 * @brief Namespace containing statistical functions
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-09-10 17:02:20 +0200 (Mon, 10 Sep 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RbMathHelper.h 1819 2012-09-10 15:02:20Z hoehna $
 */
#ifndef RbMathHelper_H
#define RbMathHelper_H

#include <vector>
#include "RandomNumberGenerator.h"
#include "RbException.h"

namespace RevBayesCore {
    
    namespace RbMath {
        
        namespace Helper {
            double          fmax2(double x, double y);
            double          fmin2(double x, double y);
        }
        
	}
}

#endif
