#include "PointMassDistribution.h"

#include <cmath>

#include "Cloneable.h"
#include "RbException.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/*
 * Dirac delta distribution Constructor
 *
 * @param v The value for the distribution
 *
 */

PointMassDistribution::PointMassDistribution(const TypedDagNode<double> *v) : ContinuousDistribution( new double( 0.0 ) ),
    val( v )
{
    
    addParameter( val );
    
    *value = val->getValue();
}


PointMassDistribution::~PointMassDistribution( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


double PointMassDistribution::cdf( void ) const
{
    
    return *value > val->getValue();
}


PointMassDistribution* PointMassDistribution::clone( void ) const
{
    
    return new PointMassDistribution( *this );
}


double PointMassDistribution::computeLnProbability( void )
{
    
    return log(*value == val->getValue());
}


void PointMassDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode* affecter)
{
    // only delegate when the toucher was our parameters
    if ( affecter == val && this->dag_node != NULL )
    {
        this->dag_node->initiateGetAffectedNodes( affected );
    }
    
}


double PointMassDistribution::getMin( void ) const
{

    return val->getValue();
}


double PointMassDistribution::getMax(void) const
{

    return val->getValue();
}


void PointMassDistribution::keepSpecialization( const DagNode* affecter )
{
    // only do this when the toucher was our parameters
    if ( affecter == val && this->dag_node != NULL )
    {
        this->dag_node->keepAffected();
    }
    
}



double PointMassDistribution::quantile(double p) const
{
    
    throw(RbException("Quantiles not defined for the Dirac delta function"));
}


void PointMassDistribution::redrawValue( void )
{
    
    *value = val->getValue();
    
}



void PointMassDistribution::restoreSpecialization( const DagNode *restorer )
{
    
    // only do this when the toucher was our parameters
    if ( restorer == val )
    {
        const double &tmp = val->getValue();
        *(this->value) = tmp;

        if ( this->dag_node != NULL )
        {
            this->dag_node->restoreAffected();
        }
        
    }
    
}


/** Swap a parameter of the distribution */
void PointMassDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == val)
    {
        val = static_cast<const TypedDagNode<double>* >( newP );
    }
    
}


void PointMassDistribution::touchSpecialization( const DagNode *toucher, bool touchAll )
{
    // only do this when the toucher was our parameters
    if ( toucher == val )
    {
        const double &tmp = val->getValue();
        *(this->value) = tmp;
        
        if ( this->dag_node != NULL )
        {
            this->dag_node->touchAffected();
        }
        
    }
    
}
