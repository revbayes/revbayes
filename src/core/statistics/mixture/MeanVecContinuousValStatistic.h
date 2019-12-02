/**
 * @file
 * This file contains the declaration of the average value statistic for a vector of doubles, 
 *
 * @brief Declaration of the MeanVecContinuousValStatistic
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: MeanVecContinuousValStatistic.h tracyh $
 */

#ifndef MeanVecContinuousValStatistic_H
#define MeanVecContinuousValStatistic_H

#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>
#include <string>

namespace RevBayesCore {
    
    class MeanVecContinuousValStatistic : public TypedFunction<double> {
        
    public:
																MeanVecContinuousValStatistic(const TypedDagNode< RbVector<double> >* v);	//!< Default constructor
																MeanVecContinuousValStatistic(const MeanVecContinuousValStatistic& t);			//!< Copy constructor
        virtual													~MeanVecContinuousValStatistic(void);											//!< Destructor
        
        MeanVecContinuousValStatistic&							operator=(const MeanVecContinuousValStatistic& t);
        
        // Basic utility functions
        MeanVecContinuousValStatistic*							clone(void) const;																//!< Clone object
        void													update(void);																	//!< Clone the function
		
    protected:
        void													swapParameterInternal(const DagNode *oldP, const DagNode *newP);				//!< Implementation of swaping parameters
        
    private:
        // members
        const TypedDagNode< RbVector<double> >*			elementVals;
    };
}


#endif

