//
//  GeneralRateGeneratorSequenceFunction.cpp
//  rb_mlandis
//
//  Created by Michael Landis on 4/3/14.
//  Copyright (c) 2014 Michael Landis. All rights reserved.
//

#include "GeneralRateGeneratorSequenceFunction.h"
#include "CharacterHistoryRateModifier.h"
#include "ConstantNode.h"
#include "RateGeneratorSequenceUsingMatrix.h"
#include "RateMatrix_JC.h"
#include "RbException.h"
#include "RbVector.h"

#include <algorithm>

using namespace RevBayesCore;

GeneralRateGeneratorSequenceFunction::GeneralRateGeneratorSequenceFunction(size_t ns, size_t nc)
    : TypedFunction<RateGeneratorSequence>( new RateGeneratorSequenceUsingMatrix(ns, nc) )
{
    homogeneousRateMatrix = new ConstantNode<RateGenerator>("homogeneousRateMatrix", new RateMatrix_JC(ns));
    rateModifiers         = NULL;
    
    this->addParameter(homogeneousRateMatrix);
    this->addParameter(rateModifiers);
    
    update();
}


GeneralRateGeneratorSequenceFunction::GeneralRateGeneratorSequenceFunction(const GeneralRateGeneratorSequenceFunction &n) : TypedFunction<RateGeneratorSequence>( n )
{
    homogeneousRateMatrix = n.homogeneousRateMatrix;
    rateModifiers         = n.rateModifiers;
}


GeneralRateGeneratorSequenceFunction::~GeneralRateGeneratorSequenceFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

GeneralRateGeneratorSequenceFunction* GeneralRateGeneratorSequenceFunction::clone( void ) const
{
    return new GeneralRateGeneratorSequenceFunction( *this );
}


void GeneralRateGeneratorSequenceFunction::update( void )
{

    const RateGenerator& rm = homogeneousRateMatrix->getValue();
    
    if (rateModifiers != NULL) {
        const RbVector<CharacterHistoryRateModifier>& rate_mods = rateModifiers->getValue();
        static_cast< RateGeneratorSequenceUsingMatrix* >(value)->setRateModifiers( &rate_mods );
    }

    // set the rate modifiers in RateGeneratorSequence!!
    static_cast< RateGeneratorSequenceUsingMatrix* >(value)->setRateMatrix( &rm );
    

    value->updateMap();
}

void GeneralRateGeneratorSequenceFunction::setRateMatrix(const TypedDagNode<RateGenerator>* r)
{
    // remove the old parameter first
    if ( homogeneousRateMatrix != NULL )
    {
        this->removeParameter( homogeneousRateMatrix );
        homogeneousRateMatrix = NULL;
    }

    // set the value
    homogeneousRateMatrix = r;

    // add the parameter
    this->addParameter( homogeneousRateMatrix );
}

void GeneralRateGeneratorSequenceFunction::setRateModifiers(const TypedDagNode<RbVector<CharacterHistoryRateModifier> >* r)
{
    // remove the old parameter first
    if ( rateModifiers != NULL )
    {
        this->removeParameter( rateModifiers );
        rateModifiers = NULL;
    }
    
    // set the value
    rateModifiers = r;
    
    // add the parameter
    this->addParameter( rateModifiers );
}


void GeneralRateGeneratorSequenceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == homogeneousRateMatrix)
    {
        homogeneousRateMatrix = static_cast<const TypedDagNode<RateGenerator>* >( newP );
    }
    else if (oldP == rateModifiers)
    {
        rateModifiers = static_cast<const TypedDagNode<RbVector<CharacterHistoryRateModifier> >* >( newP );
    }
}
