/**
 * @file
 * This file contains commonly used math functions that are used
 * in RevBayes. 
 *
 * @brief Namespace containing statistical functions
 *
 * (c) Copyright 2009- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 */
#ifndef RbMathHelper_H
#define RbMathHelper_H

namespace RevBayesCore {
    
    namespace RbMath {
        
        namespace Helper {
            double          fmax2(double x, double y);
            double          fmin2(double x, double y);
            double          reflectIntoInterval(double x, double lo, double hi);    //!< Bounce x into [lo, hi] by reflection
        }
        
	}
}

#endif
