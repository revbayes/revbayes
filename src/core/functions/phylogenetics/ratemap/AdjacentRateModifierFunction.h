//
//  AdjacentRateModifierFunction.hpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef AdjacentRateModifierFunction_h
#define AdjacentRateModifierFunction_h

#include <cstddef>
#include <iostream>

#include "CharacterHistoryRateModifier.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class AdjacentRateModifierFunction : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        AdjacentRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<std::int64_t>* w, const TypedDagNode<RbVector<RbVector<std::int64_t> > >* c, size_t ns, size_t nc);
        AdjacentRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<std::int64_t>* w, const TypedDagNode<MatrixReal>* c, size_t ns, size_t nc);
        
        AdjacentRateModifierFunction(const AdjacentRateModifierFunction& m);
        virtual ~AdjacentRateModifierFunction(void);                                                                                                  //!< Virtual destructor
        
        // public member functions
        AdjacentRateModifierFunction*                                     clone(void) const;                                                          //!< Create an independent clone
        void                                                              update(void);
        
    protected:
        void                                                              swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<double>*                                       gainFactor;
        const TypedDagNode<double>*                                       lossFactor;
        const TypedDagNode<std::int64_t>*                                         width;
        const TypedDagNode<RbVector<RbVector<std::int64_t> > >*                   context_array;
        const TypedDagNode<MatrixReal>*                                   context_matrix;
        
        std::string context_type;
        
        
    };
}

#endif /* AdjacentRateModifierFunction_hpp */
