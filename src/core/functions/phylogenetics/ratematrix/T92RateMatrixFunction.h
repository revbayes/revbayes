/**
 * @file
 * This file contains the declaration of the GTR rate matrix function class.
 * This class is derived from the function class and is used to
 * compute the rate matrix of a general time reversible (GTR) Markov chain.
 *
 * @brief Declaration of the GTR rate matrix function.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface Function
 *
 * $Id$
 */



#ifndef T92RateMatrixFunction_H
#define T92RateMatrixFunction_H

#include "RateMatrix_GTR.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class T92RateMatrixFunction : public TypedFunction<RateGenerator> {
        
    public:
        T92RateMatrixFunction(const TypedDagNode<double > *eqGc, const TypedDagNode< double > *tstv);
        virtual                                            ~T92RateMatrixFunction(void);                                                            //!< Virtual destructor
        
        // public member functions
        T92RateMatrixFunction*                              clone(void) const;                                                                      //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< double >*                       equilibriumGc;
        const TypedDagNode< double >*                       transitionTransversionRate;
        
    };
    
}

#endif
