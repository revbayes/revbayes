/**
 * @file
 * This file contains the declaration of the PoMo root frequencies funcion class.
 * This class is derived from the function class and is used to
 * compute the root frequencies of the PoMo model.
 *
 * @brief Declaration of the PoMo root frequencies function.
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



#ifndef PoMoRootFrequenciesFunction_H
#define PoMoRootFrequenciesFunction_H

#include <vector>
#include <cstdint>

#include "Simplex.h"
#include "TypedFunction.h"
#include "RateGenerator.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class PoMoRootFrequenciesFunction : public TypedFunction< Simplex > {
        
    public:
        PoMoRootFrequenciesFunction(const TypedDagNode< Simplex > *fnrf, const TypedDagNode< double > *fopar, const TypedDagNode< RbVector<double> > *mr, const TypedDagNode< std::int64_t > *ps);
        
        PoMoRootFrequenciesFunction(const TypedDagNode< Simplex > *fnrf, const TypedDagNode< double > *fopar, const TypedDagNode< RateGenerator > *mm, const TypedDagNode< std::int64_t > *ps);

        virtual                                            ~PoMoRootFrequenciesFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        PoMoRootFrequenciesFunction*                        clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< Simplex >*                      fixedNucleotideRootFrequencies;
        const TypedDagNode< double >*                       frequencyOfPolymorphismsAtTheRoot;
        const TypedDagNode< RbVector<double> >*             mutationRates;
        const TypedDagNode< RateGenerator >*                mutationMatrix;
        const TypedDagNode< std::int64_t >*                          populationSize;
        bool                                                useMutationMatrix;

        // Private member functions
        std::vector<double> setMutationRates(const RateGenerator& mm) ;

        
    };
    
}

#endif
