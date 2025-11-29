/**
 * @file
 * This file contains the implementation of the ChromosomesPloidyRateMatrixFunction class.
 * This class is derived from the function class and is used to
 * create the chromosomes rate matrix.
 *
 * @brief Implementation of the ChromosomesPloidyRateMatrixFunction.
 *
 * (c) Copyright 2014- under GPL version 3
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 */

#include "ChromosomesPloidyRateMatrixFunction.h"

#include "RateMatrix_ChromosomesPloidy.h"
#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

ChromosomesPloidyRateMatrixFunction::ChromosomesPloidyRateMatrixFunction(const TypedDagNode<std::int64_t> *n, const TypedDagNode<double> *g_d, const TypedDagNode<double> *d_d, const TypedDagNode<double> *r_d, const TypedDagNode<double> *e_d, const TypedDagNode<double> *g_dl, const TypedDagNode<double> *d_dl, const TypedDagNode<double> *g_p, const TypedDagNode<double> *d_p, const TypedDagNode<double> *r_p, const TypedDagNode<double> *e_p, const TypedDagNode<double> *g_pl, const TypedDagNode<double> *d_pl, const TypedDagNode<double> *b) : TypedFunction<RateGenerator>( new RateMatrix_ChromosomesPloidy(n->getValue()) ), gamma_d( g_d ), delta_d( d_d ), rho_d( r_d ), eta_d( e_d ), gamma_dl( g_dl ), delta_dl( d_dl ), gamma_p( g_p ), delta_p( d_p ), rho_p( r_p ), eta_p( e_p ), gamma_pl( g_pl ), delta_pl( d_pl ), beta( b ) {

    addParameter( gamma_d );
    addParameter( delta_d );
    addParameter( rho_d );
    addParameter( eta_d );
    addParameter( gamma_dl );
    addParameter( delta_dl );
    addParameter( gamma_p );
    addParameter( delta_p );
    addParameter( rho_p );
    addParameter( eta_p );
    addParameter( gamma_pl );
    addParameter( delta_pl );
    addParameter( beta );
    
    update();
}



ChromosomesPloidyRateMatrixFunction::~ChromosomesPloidyRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



ChromosomesPloidyRateMatrixFunction* ChromosomesPloidyRateMatrixFunction::clone( void ) const {
    return new ChromosomesPloidyRateMatrixFunction( *this );
}


void ChromosomesPloidyRateMatrixFunction::update( void ) {
    double ga_d = gamma_d->getValue();
    double de_d = delta_d->getValue();
    double rh_d = rho_d->getValue();
    double ga_dl = gamma_dl->getValue();
    double de_dl = delta_dl->getValue();
    double e_d = eta_d->getValue();
    double ga_p = gamma_p->getValue();
    double de_p = delta_p->getValue();
    double rh_p = rho_p->getValue();
    double ga_pl = gamma_pl->getValue();
    double de_pl = delta_pl->getValue();
    double e_p = eta_p->getValue();
    double be = beta->getValue();

    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setGamma_d( ga_d );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setDelta_d( de_d );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setRho_d( rh_d );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setEta_d( e_d );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setGamma_dl( ga_dl );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setDelta_dl( de_dl );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setGamma_p( ga_p );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setDelta_p( de_p );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setRho_p( rh_p );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setEta_p( e_p );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setGamma_pl( ga_pl );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setDelta_pl( de_pl );
    static_cast< RateMatrix_ChromosomesPloidy* >(value)->setBeta( be );

    static_cast< RateMatrix_ChromosomesPloidy* >(value)->update();
    
}



void ChromosomesPloidyRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == gamma_d) 
    {
        gamma_d = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == delta_d) 
    {
        delta_d = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho_d) 
    {
        rho_d = static_cast<const TypedDagNode<double>* >( newP );
    } else if (oldP == eta_d) 
    {
        eta_d = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == delta_dl) 
    {
        delta_dl = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == gamma_dl) 
    {
        gamma_dl = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == gamma_p) 
    {
        gamma_p = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == delta_p) 
    {
        delta_p = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == rho_p) 
    {
        rho_p = static_cast<const TypedDagNode<double>* >( newP );
    } else if (oldP == eta_p) 
    {
        eta_p = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == delta_pl) 
    {
        delta_pl = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == gamma_pl) 
    {
        gamma_pl = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == beta) 
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}



