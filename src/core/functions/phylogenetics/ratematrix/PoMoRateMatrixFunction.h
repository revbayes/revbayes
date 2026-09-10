/**
 * @file
 * This file contains the declaration of the PoMo rate matrix function class.
 * This class is derived from the function class and is used to
 * compute the rate matrix of a general time reversible (PoMo) Markov chain.
 *
 * @brief Declaration of the PoMo rate matrix function.
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



#ifndef PoMoRateMatrixFunction_H
#define PoMoRateMatrixFunction_H

#include <vector>

#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class PoMoRateMatrixFunction : public TypedFunction<RateGenerator> {

    public:
        PoMoRateMatrixFunction(const TypedDagNode< std::int64_t > *ps, const TypedDagNode< RbVector<double> > *mr, const TypedDagNode< RbVector<double>  > *sc);
        PoMoRateMatrixFunction(const TypedDagNode< std::int64_t > *ps, const TypedDagNode< RateGenerator > *mm, const TypedDagNode< RbVector<double>  > *sc);

        virtual                                            ~PoMoRateMatrixFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoRateMatrixFunction*                             clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members

        const TypedDagNode< std::int64_t >*                          populationSize;
        const TypedDagNode< RbVector<double> >*             mutationRates;
        const TypedDagNode< RateGenerator >*                mutationMatrix;
        const TypedDagNode< RbVector<double> >*             selectionCoefficients;
        bool                                                useMutationMatrix;
        
        // private member functions
        std::vector<double> setMutationRates(const RateGenerator& mm);

    };
    
}

#endif
