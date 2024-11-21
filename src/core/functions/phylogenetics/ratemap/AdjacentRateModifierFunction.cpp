//
//  AdjacentRateModifierFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 2/3/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#include "AdjacentRateModifierFunction.h"

#include <string>
#include <cstdint>

#include "AdjacentRateModifier.h"
#include "TypedDagNode.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

AdjacentRateModifierFunction::AdjacentRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<std::int64_t>* w, const TypedDagNode<RbVector<RbVector<std::int64_t> > >* c, size_t ns, size_t nc) : TypedFunction<CharacterHistoryRateModifier>( new AdjacentRateModifier(ns, nc) ),
    gainFactor(gf),
    lossFactor(lf),
    width(w),
    context_matrix(NULL),
    context_array(c),
    context_type("width")
{
    if (context_array != NULL) {
        context_type = "array";
    }
    
    // add the parameters as parents
    addParameter(gainFactor);
    addParameter(lossFactor);
    addParameter(width);
    addParameter(context_array);
    
    update();
}

AdjacentRateModifierFunction::AdjacentRateModifierFunction(const TypedDagNode<double>* gf, const TypedDagNode<double>* lf, const TypedDagNode<std::int64_t>* w, const TypedDagNode<MatrixReal>* c, size_t ns, size_t nc) : TypedFunction<CharacterHistoryRateModifier>( new AdjacentRateModifier(ns, nc) ),
gainFactor(gf),
lossFactor(lf),
width(w),
context_array(NULL),
context_matrix(c),
context_type("width")
{
    if (context_matrix != NULL) {
        context_type = "matrix";
    }
    
    // add the parameters as parents
    addParameter(gainFactor);
    addParameter(lossFactor);
    addParameter(width);
    addParameter(context_matrix);
    
    update();
}

AdjacentRateModifierFunction::AdjacentRateModifierFunction(const AdjacentRateModifierFunction& m) : TypedFunction<CharacterHistoryRateModifier>( m )
{
    gainFactor = m.gainFactor;
    lossFactor = m.lossFactor;
    width = m.width;
    context_matrix = m.context_matrix;
    context_array = m.context_array;
    context_type = m.context_type;
}


AdjacentRateModifierFunction::~AdjacentRateModifierFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}





AdjacentRateModifierFunction* AdjacentRateModifierFunction::clone( void ) const
{
    return new AdjacentRateModifierFunction( *this );
}


void AdjacentRateModifierFunction::update( void )
{
    
    double gf = gainFactor->getValue();
    static_cast<AdjacentRateModifier*>(value)->setGainFactor(gf);
    
    double lf = lossFactor->getValue();
    static_cast<AdjacentRateModifier*>(value)->setLossFactor(lf);
    
    if (context_type == "width")
    {
        size_t w = width->getValue();
        static_cast<AdjacentRateModifier*>(value)->setWidth(w);
    }
    else if (context_type=="array")
    {
        RbVector<RbVector<std::int64_t> > c = context_array->getValue();
        static_cast<AdjacentRateModifier*>(value)->setContextMatrix(c);
    }
    else if (context_type=="matrix")
    {
        MatrixReal c = context_matrix->getValue();
        static_cast<AdjacentRateModifier*>(value)->setContextMatrix(c);
    }
    
}



void AdjacentRateModifierFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == gainFactor)
    {
        gainFactor = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == lossFactor)
    {
        lossFactor = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == width)
    {
        width = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    else if (oldP == context_array)
    {
        context_array = static_cast<const TypedDagNode<RbVector<RbVector<std::int64_t> > >* >( newP );
    }
    else if (oldP == context_matrix)
    {
        context_matrix = static_cast<const TypedDagNode<MatrixReal>* >( newP );
    }
}
