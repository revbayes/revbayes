/**
 * @file
 * This file contains the implementation of the ChromosomesRateMatrixFunction class.
 * This class is derived from the function class and is used to
 * create the chromosomes rate matrix.
 *
 * @brief Implementation of the ChromosomesRateMatrixFunction.
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 */

#include "ChromosomesRateMatrixFunction.h"

#include "RateMatrix_Chromosomes.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

ChromosomesRateMatrixFunction::ChromosomesRateMatrixFunction(size_t n, const TypedDagNode<double> *g, const TypedDagNode<double> *d, const TypedDagNode<double> *r, const TypedDagNode<double> *e, const TypedDagNode<double> *g_l, const TypedDagNode<double> *d_l, AbstractRateMatrix::METHOD m) : TypedFunction<RateGenerator>( new RateMatrix_Chromosomes(n, m) ),
    gamma( g ),
    delta( d ),
    rho( r ),
    eta( e ),
    gamma_l( g_l ),
    delta_l( d_l ) 
{

    addParameter( gamma );
    addParameter( delta );
    addParameter( rho );
    addParameter( eta );
    addParameter( gamma_l );
    addParameter( delta_l );
    
    update();
}



ChromosomesRateMatrixFunction* ChromosomesRateMatrixFunction::clone( void ) const 
{
    return new ChromosomesRateMatrixFunction( *this );
}


void ChromosomesRateMatrixFunction::update( void ) 
{
    double ga = gamma->getValue();
    double de = delta->getValue();
    double rh = rho->getValue();
    double ga_l = gamma_l->getValue();
    double de_l = delta_l->getValue();
    double e = eta->getValue();

    static_cast< RateMatrix_Chromosomes* >(value)->setGamma( ga );
    static_cast< RateMatrix_Chromosomes* >(value)->setDelta( de );
    static_cast< RateMatrix_Chromosomes* >(value)->setRho( rh );
    static_cast< RateMatrix_Chromosomes* >(value)->setEta( e );
    static_cast< RateMatrix_Chromosomes* >(value)->setGamma_l( ga_l );
    static_cast< RateMatrix_Chromosomes* >(value)->setDelta_l( de_l );

    static_cast< RateMatrix_Chromosomes* >(value)->update();
    
}



void ChromosomesRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) 
{
    
    if (oldP == gamma) 
    {
        gamma = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == delta)
    {
        delta = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode<double>* >( newP );
    } 
    
    if (oldP == eta)
    {
        eta = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == delta_l)
    {
        delta_l = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == gamma_l)
    {
        gamma_l = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}



