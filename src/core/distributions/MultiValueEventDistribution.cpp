#include "MultiValueEventDistribution.h"

#include <cstddef>
#include <string>

#include "DagNode.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbVector.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

MultiValueEventDistribution::MultiValueEventDistribution(TypedDistribution<long> *ep, const std::vector< TypedDistribution<double> *> &vp, const std::vector< std::string > &n, const std::vector< long > &min) : TypedDistribution< MultiValueEvent >( new MultiValueEvent() ),
    event_prior( ep ),
    min_events( min ),
    names( n ),
    value_priors( vp )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& even_prior_pars = event_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = even_prior_pars.begin(); it != even_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    // add the parameters of the distribution
    for ( size_t i=0; i<value_priors.size(); ++i )
    {
        const std::vector<const DagNode*>& value_prior_pars = value_priors[i]->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    }
    
    simulate();
}



MultiValueEventDistribution::MultiValueEventDistribution(const MultiValueEventDistribution &d) : TypedDistribution<MultiValueEvent>( d ),
    event_prior( d.event_prior->clone() ),
    min_events( d.min_events ),
    names( d.names ),
    value_priors()
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    const std::vector<const DagNode*>& even_prior_pars = event_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = even_prior_pars.begin(); it != even_prior_pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    // add the parameters of the distribution
    for ( size_t i=0; i<d.value_priors.size(); ++i )
    {
        // first we need to clone the base distribution
        value_priors.push_back( d.value_priors[i]->clone() );
        
        const std::vector<const DagNode*>& value_prior_pars = value_priors[i]->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    }
    
}



MultiValueEventDistribution& MultiValueEventDistribution::operator=(const MultiValueEventDistribution &d)
{
    
    if ( this != &d )
    {
        TypedDistribution<MultiValueEvent>::operator=( d );
        
        delete event_prior;
        for ( size_t i=0; i<value_priors.size(); ++i )
        {
            TypedDistribution<double> *tmp_dist = value_priors[i];
            delete tmp_dist;
        }
        
        event_prior         = d.event_prior->clone();
        min_events          = d.min_events;
        names               = d.names;
        value_priors.clear();
        
        // add the parameters of the distribution
        for ( size_t i=0; i<d.value_priors.size(); ++i )
        {
            // first we need to clone the base distribution
            value_priors.push_back( d.value_priors[i]->clone() );
        }
        
    }
    
    return *this;
}


MultiValueEventDistribution::~MultiValueEventDistribution( void )
{
    
    delete event_prior;
    for ( size_t i=0; i<value_priors.size(); ++i )
    {
        TypedDistribution<double> *tmp_dist = value_priors[i];
        delete tmp_dist;
    }
    
}




MultiValueEventDistribution* MultiValueEventDistribution::clone( void ) const
{
    
    return new MultiValueEventDistribution( *this );
}



double MultiValueEventDistribution::computeLnProbability( void )
{
    event_prior->setValue( new long( value->getNumberOfEvents() ) );
    
    double ln_prob = event_prior->computeLnProbability();
    
    for (int j = 0; j < value_priors.size(); ++j)
    {

        const std::vector<double> &these_values = this->value->getValues(j);
        
        long this_num_values = these_values.size();
        
        if ( (value->getNumberOfEvents()+min_events[j]) != this_num_values )
        {
            return RbConstants::Double::neginf;
        }

        for (int i = 0; i < this_num_values; ++i)
        {

            // we also need to multiply with the probability of the value for this table
            value_priors[j]->setValue( new double( these_values[i] ) );
            ln_prob += value_priors[j]->computeLnProbability();
            
        }
        
    }
    
    return ln_prob;
}



void MultiValueEventDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    
    if ( n == "getValues" )
    {
        const std::string &val_name = static_cast<const TypedDagNode<std::string> * >( args[0] )->getValue();
        rv = value->getValues( val_name );
        
    }
    else if ( n == "getRealValues" )
    {
        const std::string &val_name = static_cast<const TypedDagNode<std::string> * >( args[0] )->getValue();
        rv = value->getValues( val_name );
    }
    else if ( n == "getRealPosValues" )
    {
        const std::string &val_name = static_cast<const TypedDagNode<std::string> * >( args[0] )->getValue();
        rv = value->getValues( val_name );
    }
    else
    {
        throw RbException("The multi-value event does not have a member method called '" + n + "'.");
    }
    
}



void MultiValueEventDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, long &rv) const
{
    
    if ( n == "getNumberOfEvents" )
    {
        rv = value->getNumberOfEvents();
        
    }
    //    else if ( n == "getNumberOfEvents" )
    //    {
    //
    //
    //    }
    else
    {
        throw RbException("The multi-value event does not have a member method called '" + n + "'.");
    }
    
}

const std::vector<long>& MultiValueEventDistribution::getMinimumNumberOfEvents(void) const
{
    return min_events;
}


const std::vector< TypedDistribution<double>* >& MultiValueEventDistribution::getValuePriors(void) const
{
    
    return value_priors;
}



void MultiValueEventDistribution::simulate()
{
    // clear the current value
    this->value->clear();
    
    // draw a number of events
    event_prior->redrawValue();
    value->setNumberOfEvents( event_prior->getValue() );
    
    for (int j = 0; j < value_priors.size(); ++j)
    {
        
        long this_num_values = min_events[j] + value->getNumberOfEvents();
        
        std::vector<double>      these_values   = std::vector<double>(this_num_values, 0.0);
        const std::string       &this_name      = names[j];
        TypedDistribution<double>   *this_prior = value_priors[j];

        for (long i = 0; i < this_num_values; ++i)
        {
        
            // we also need to multiply with the probability of the value for this table
            this_prior->redrawValue();
            these_values[i] = this_prior->getValue();
            
        }
        
        this->value->addValues(these_values, this_name);
        
    }
    
    
}


void MultiValueEventDistribution::redrawValue( void )
{
    
    simulate();
    
}



/** Swap a parameter of the distribution */
void MultiValueEventDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    bool found = false;
    try {
        event_prior->swapParameter(oldP,newP);
        // if the statement succeeded and didn't throw an error, then the distribution had this parameter
        found = true;
    }
    catch (RbException &e)
    {
        // do nothing because we actually do not know who had the parameter
    }
    
    for (int j = 0; j < value_priors.size(); ++j)
    {
        
        TypedDistribution<double> *this_prior = value_priors[j];
        try {
            this_prior->swapParameter(oldP,newP);
            // if the statement succeeded and didn't throw an error, then the distribution had this parameter
            found = true;
        }
        catch (RbException &e)
        {
            // do nothing because we actually do not know who had the parameter
        }
    
    }
    
    
    if ( found == false )
    {
        throw RbException("Could not find the distribution parameter to be swapped: " + oldP->getName() + " to " + newP->getName()) ;
    }
    
}
