#include "AutocorrelatedEventDistribution.h"

#include <cstddef>
#include <cstdint>
#include <string>

#include "DagNode.h"
#include "DistributionNormal.h"
#include "DistributionLognormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbVector.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

AutocorrelatedEventDistribution::AutocorrelatedEventDistribution(TypedDistribution<std::int64_t> *ep, const std::vector< TypedDistribution<double> *>& vp, const std::vector< Autocorrelation >& ac, const std::vector< std::string >& ac_dep_var, const TypedDagNode< RbVector<double> >* ac_sd, const std::vector< std::string >& n, const std::vector< std::int64_t >& min, const std::string& sort_var) : TypedDistribution< MultiValueEvent >( new MultiValueEvent() ),
    event_prior( ep ),
    min_events( min ),
    names( n ),
    value_priors( vp ),
    autocorrelation_types( ac ),
    name_of_var_to_sort_by( sort_var ),
    autocorrelation_sigmas( ac_sd )
{

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& even_prior_pars = event_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = even_prior_pars.begin(); it != even_prior_pars.end(); ++it)
    {
        addParameter( *it );
    }
    
    // add the parameters of the distribution
    for ( size_t i=0; i<value_priors.size(); ++i )
    {
        const std::vector<const DagNode*>& value_prior_pars = value_priors[i]->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
        {
            addParameter( *it );
        }
    }
    
    // add the sd parameters of the autocorrelated distribution
    addParameter( autocorrelation_sigmas);
    
    // initialize the index variables
    index_of_var_to_sort_by = -1;
    for ( int i=0; i<value_priors.size(); ++i )
    {
        if (names[i] == name_of_var_to_sort_by)
        {
            index_of_var_to_sort_by = i;
            break;
        }
    }
    autocorrelation_time_indeces = std::vector< int >(value_priors.size(), -1);
    for ( int i=0; i<value_priors.size(); ++i )
    {
        if ( autocorrelation_types[i] != NONE )
        {
            for ( int j=0; j<value_priors.size(); ++j )
            {
                if (names[j] == ac_dep_var[i])
                {
                    autocorrelation_time_indeces[i] = j;
                    break;
                }
            }
        }
    }
    
    simulate();
}



AutocorrelatedEventDistribution::AutocorrelatedEventDistribution(const AutocorrelatedEventDistribution &d) : TypedDistribution<MultiValueEvent>( d ),
    event_prior( d.event_prior->clone() ),
    min_events( d.min_events ),
    names( d.names ),
    value_priors(),
    autocorrelation_types( d.autocorrelation_types ),
    name_of_var_to_sort_by( d.name_of_var_to_sort_by ),
    index_of_var_to_sort_by( d.index_of_var_to_sort_by ),
    autocorrelation_sigmas( d.autocorrelation_sigmas ),
    autocorrelation_time_indeces( d.autocorrelation_time_indeces )
{
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    const std::vector<const DagNode*>& even_prior_pars = event_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = even_prior_pars.begin(); it != even_prior_pars.end(); ++it)
    {
        addParameter( *it );
    }
    
    // add the parameters of the distribution
    for ( size_t i=0; i<d.value_priors.size(); ++i )
    {
        // first we need to clone the base distribution
        value_priors.push_back( d.value_priors[i]->clone() );
        
        const std::vector<const DagNode*>& value_prior_pars = value_priors[i]->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = value_prior_pars.begin(); it != value_prior_pars.end(); ++it)
        {
            addParameter( *it );
        }
    }
    
    // add the sd parameters of the autocorrelated distribution
    addParameter( autocorrelation_sigmas );
    
}



AutocorrelatedEventDistribution& AutocorrelatedEventDistribution::operator=(const AutocorrelatedEventDistribution &d)
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
        autocorrelation_types           = d.autocorrelation_types;
        name_of_var_to_sort_by          = d.name_of_var_to_sort_by;
        index_of_var_to_sort_by         = d.index_of_var_to_sort_by;
        autocorrelation_sigmas          = d.autocorrelation_sigmas;
        autocorrelation_time_indeces    = d.autocorrelation_time_indeces;
        
        // add the parameters of the distribution
        for ( size_t i=0; i<d.value_priors.size(); ++i )
        {
            // first we need to clone the base distribution
            value_priors.push_back( d.value_priors[i]->clone() );
        }
        
    }
    
    return *this;
}


AutocorrelatedEventDistribution::~AutocorrelatedEventDistribution( void )
{
    
    delete event_prior;
    for ( size_t i=0; i<value_priors.size(); ++i )
    {
        TypedDistribution<double> *tmp_dist = value_priors[i];
        delete tmp_dist;
    }
    
}




AutocorrelatedEventDistribution* AutocorrelatedEventDistribution::clone( void ) const
{
    
    return new AutocorrelatedEventDistribution( *this );
}



double AutocorrelatedEventDistribution::computeLnProbability( void )
{
    
    // first, we check the sorting (-1 if we don't sort)
    if ( index_of_var_to_sort_by >= 0 )
    {
        const std::vector<double> &these_values = this->value->getValues( index_of_var_to_sort_by );
            
        std::int64_t this_num_values = these_values.size();
        
        for (int i = 1; i < this_num_values; ++i)
        {

            // we also need to multiply with the probability of the value for this table
            if ( these_values[i-1] > these_values[i] )
            {
                return RbConstants::Double::neginf;
            }
            
        }
    }
    
    
    // compute the prior probability for the number of events
    event_prior->setValue( new std::int64_t( value->getNumberOfEvents() ) );
    double ln_prob = event_prior->computeLnProbability();
    
    
    // compute the priors for each new valuea at each event
    for (int j = 0; j < value_priors.size(); ++j)
    {

        const std::vector<double> &these_values = this->value->getValues(j);
        
        std::int64_t this_num_values = these_values.size();
        
        if ( (value->getNumberOfEvents()+min_events[j]) != this_num_values )
        {
            return RbConstants::Double::neginf;
        }
        
        const std::vector<double> &time_values = this->value->getValues( autocorrelation_time_indeces[j] );

        for (int i = 0; i < this_num_values; ++i)
        {

            // we also need to multiply with the probability of the value for this table
            // first, if this variable is uncorrelated or the first value
            if ( i == 0 || autocorrelation_types[j] == NONE )
            {
                value_priors[j]->setValue( new double( these_values[i] ) );
                ln_prob += value_priors[j]->computeLnProbability();
            }
            // if this variable is correlate via a normal distribution
            else if ( autocorrelation_types[j] == ACN )
            {
                std::int64_t time_index = i - min_events[j] + min_events[autocorrelation_time_indeces[j]];

                double mean = these_values[i-1];
                double val  = these_values[i];
                double dt   = time_values[time_index] - (time_index > 0 ? time_values[time_index-1] : 0.0);
                double sd   = autocorrelation_sigmas->getValue()[j] * dt ;
                ln_prob += RbStatistics::Normal::lnPdf(mean, sd, val);
                
            }
            // if this variable is correlated via a lognormal distribution
            else if ( autocorrelation_types[j] == ACLN )
            {
                std::int64_t time_index = i - min_events[j] + min_events[autocorrelation_time_indeces[j]];
                
                double mean = log( these_values[i-1] );
                double val  = these_values[i];
                double dt   = time_values[time_index] - (time_index > 0 ? time_values[time_index-1] : 0.0);
                double sd   = autocorrelation_sigmas->getValue()[j] * dt ;
                
                ln_prob += RbStatistics::Lognormal::lnPdf(mean, sd, val);
                
            }
            
        }
        
    }
    
    return ln_prob;
}



void AutocorrelatedEventDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
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



void AutocorrelatedEventDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, std::int64_t &rv) const
{
    
    if ( n == "getNumberOfEvents" )
    {
        rv = value->getNumberOfEvents();
        
    }
    else
    {
        throw RbException("The multi-value event does not have a member method called '" + n + "'.");
    }
    
}

const std::vector<std::int64_t>& AutocorrelatedEventDistribution::getMinimumNumberOfEvents(void) const
{
    return min_events;
}


const std::vector< TypedDistribution<double>* >& AutocorrelatedEventDistribution::getValuePriors(void) const
{
    
    return value_priors;
}



bool AutocorrelatedEventDistribution::isAutocorrelated(size_t i) const
{
    return autocorrelation_types[i] != NONE;
}



bool AutocorrelatedEventDistribution::isSorted(size_t i) const
{
    return i == index_of_var_to_sort_by;
}
 


void AutocorrelatedEventDistribution::simulate()
{
    // clear the current value
    this->value->clear();
    
    // draw a number of events
    event_prior->redrawValue();
    value->setNumberOfEvents( event_prior->getValue() );
    
    // we need to specify an ordering for the simulation
    // it is important that we simulate the times first if there are autocorrelated variables
    std::vector<size_t> simulation_ordering;
    std::vector<bool>   used = std::vector<bool>(value_priors.size(), false);
    // first, find and add the times
    for (size_t j = 0; j < value_priors.size(); ++j)
    {
        if ( names[j] == "time" || names[j] == "times" )
        {
            simulation_ordering.push_back(j);
            used[j] = true;
        }
    }
    // second, add the rest
    for (size_t j = 0; j < value_priors.size(); ++j)
    {
        if ( used[j] == false )
        {
            simulation_ordering.push_back(j);
        }
    }

    for (size_t k = 0; k < value_priors.size(); ++k)
    {
        size_t j = simulation_ordering[k];
        
        std::int64_t this_num_values = min_events[j] + value->getNumberOfEvents();
        
        std::vector<double>         these_values    = std::vector<double>(this_num_values, 0.0);
        const std::string&          this_name       = names[j];
        TypedDistribution<double>*  this_prior      = value_priors[j];
        
        const std::vector<double> &time_values      = this->value->getValues( autocorrelation_time_indeces[j] );


        for (int i = 0; i < this_num_values; ++i)
        {
        
            // first, if this variable is uncorrelated or the first value
            if ( i == 0 || autocorrelation_types[j] == NONE )
            {
                this_prior->redrawValue();
                these_values[i] = this_prior->getValue();
            }
            // if this variable is correlate via a normal distribution
            else if ( autocorrelation_types[j] == ACN )
            {
                std::int64_t time_index = i - min_events[j] + min_events[autocorrelation_time_indeces[j]];

                
                double mean     = these_values[i-1];
                double dt       = time_values[time_index] - (time_index > 0 ? time_values[time_index-1] : 0.0);
                double sd       = autocorrelation_sigmas->getValue()[j] * dt ;
                these_values[i] = RbStatistics::Normal::rv(mean, sd, *GLOBAL_RNG);
                
            }
            // if this variable is correlated via a lognormal distribution
            else if ( autocorrelation_types[j] == ACLN )
            {
                std::int64_t time_index = i - min_events[j] + min_events[autocorrelation_time_indeces[j]];
                
                double mean     = log( these_values[i-1] );
                double dt       = time_values[time_index] - (time_index > 0 ? time_values[time_index-1] : 0.0);
                double sd       = autocorrelation_sigmas->getValue()[j] * dt ;
                these_values[i] = exp(RbStatistics::Normal::rv(mean, sd, *GLOBAL_RNG));
                
            }
            
        }
        
        if ( name_of_var_to_sort_by == this_name )
        {
            std::sort(these_values.begin(), these_values.end());
        }
        
        this->value->addValues(these_values, this_name);
        
    }
    
    
}


void AutocorrelatedEventDistribution::redrawValue( void )
{
    
    simulate();
    
}



/** Swap a parameter of the distribution */
void AutocorrelatedEventDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
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
    
    if ( oldP == autocorrelation_sigmas )
    {
        autocorrelation_sigmas = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
        found = true;
    }
    
    
    if ( found == false )
    {
        throw RbException("Could not find the distribution parameter to be swapped: " + oldP->getName() + " to " + newP->getName()) ;
    }
    
}
