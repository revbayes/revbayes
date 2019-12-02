//
//  AdjacentRateModifierFunction.hpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef AdjacentRateModifierFunction_h
#define AdjacentRateModifierFunction_h

#include "CharacterHistoryRateModifier.h"
#include "MatrixReal.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <iostream>
#include <vector>

namespace RevBayesCore {
    
    class AdjacentRateModifierFunction : public TypedFunction<CharacterHistoryRateModifier> {
        
    public:
        AdjacentRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<long>* w, const TypedDagNode<RbVector<RbVector<long> > >* c, size_t ns, size_t nc);
        AdjacentRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<long>* w, const TypedDagNode<MatrixReal>* c, size_t ns, size_t nc);
        
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
        const TypedDagNode<long>*                                         width;
        const TypedDagNode<RbVector<RbVector<long> > >*                   context_array;
        const TypedDagNode<MatrixReal>*                                   context_matrix;
        
        std::string context_type;
        
        
    };
}

#endif /* AdjacentRateModifierFunction_hpp */
