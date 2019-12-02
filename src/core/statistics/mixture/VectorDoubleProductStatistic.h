/**
 * @file
 * This file contains the declaration of the average value statistic for a vector of doubles, 
 *
 * @brief Declaration of the VectorDoubleProductStatistic
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: VectorDoubleProductStatistic.h tracyh $
 */

#ifndef VectorDoubleProductStatistic_H
#define VectorDoubleProductStatistic_H

#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>
#include <string>

namespace RevBayesCore {
    
    class VectorDoubleProductStatistic : public TypedFunction< RbVector<double> > {
        
    public:
		VectorDoubleProductStatistic(const TypedDagNode< RbVector<double> >* v, const TypedDagNode< double >* r);	//!< Default constructor
		VectorDoubleProductStatistic(const VectorDoubleProductStatistic& t);			//!< Copy constructor
        virtual													~VectorDoubleProductStatistic(void);											//!< Destructor
        
        VectorDoubleProductStatistic&							operator=(const VectorDoubleProductStatistic& t);
        
        // Basic utility functions
        VectorDoubleProductStatistic*							clone(void) const;																//!< Clone object
        void													update(void);																	//!< Clone the function
		
    protected:
        void													swapParameterInternal(const DagNode *oldP, const DagNode *newP);				//!< Implementation of swaping parameters
        
    private:
        // members
        const TypedDagNode< RbVector<double> >*			elementVals;
		const TypedDagNode<double>*							multiplier;
    };
}


#endif

