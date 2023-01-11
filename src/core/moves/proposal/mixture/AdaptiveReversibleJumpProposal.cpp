#include "AdaptiveReversibleJumpProposal.h"
#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "ReversibleJumpMixtureConstantDistribution.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
AdaptiveReversibleJumpProposal::AdaptiveReversibleJumpProposal( StochasticNode<double> *n, size_t n0, size_t c0, size_t ue ) : Proposal(),
    variable( n ),
    stored_value( 0 ),
    stored_index( 0 ),
    wait_before_learning( n0 ),
    wait_before_using ( c0 ),
    updates_every( ue ),
    num_tried ( 0 ),
    proposal_distribution( NORMAL )
{
    if (wait_before_using < wait_before_learning)
    {
        throw RbException("Cannot delay learning empirical proposal distribution longer than using empirical proposal distribution in adaptive RJ-switch move");
    }
    
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
AdaptiveReversibleJumpProposal::AdaptiveReversibleJumpProposal( const AdaptiveReversibleJumpProposal &p ) : Proposal(),
    variable( p.variable ),
    stored_index( p.stored_index ),
    stored_value( p.stored_value ),
    wait_before_learning( p.wait_before_learning ),
    wait_before_using ( p.wait_before_using ),
    updates_every ( p.updates_every ),
    num_tried ( p.num_tried )
{
    // tell the base class to add the node
    addNode( variable );
    
}


/**
 * Destructor
 *
 */
AdaptiveReversibleJumpProposal::~AdaptiveReversibleJumpProposal( )
{
    
}

/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void RevBayesCore::AdaptiveReversibleJumpProposal::cleanProposal( void )
{
    
}


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
RevBayesCore::AdaptiveReversibleJumpProposal& RevBayesCore::AdaptiveReversibleJumpProposal::operator=( const AdaptiveReversibleJumpProposal &p )
{
    
    if ( this != &p )
    {
        Proposal::Cloneable::operator=( p );
    
        
        stored_value    = p.stored_value;
        variable        = p.variable;
        stored_index    = p.stored_index;
        
    }
    
    return *this;
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
RevBayesCore::AdaptiveReversibleJumpProposal* RevBayesCore::AdaptiveReversibleJumpProposal::clone( void ) const
{
    
    return new AdaptiveReversibleJumpProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& RevBayesCore::AdaptiveReversibleJumpProposal::getProposalName( void ) const
{
    static std::string name = "Adaptive RJ-Switch";
    
    return name;
}


double RevBayesCore::AdaptiveReversibleJumpProposal::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * The reversible jump proposal switches the current "dimension".
 *
 * \return The hastings ratio.
 */
double RevBayesCore::AdaptiveReversibleJumpProposal::doProposal( void )
{
    
    
    double& v = variable->getValue();
    ReversibleJumpMixtureConstantDistribution<double> &d = static_cast< ReversibleJumpMixtureConstantDistribution<double>& >( variable->getDistribution() );
    
    // copy value
    stored_value = v;
    stored_index = d.getCurrentIndex();
    
    ++num_tried;
    
    if ( num_tried == 1 )
    // First time using move, setting up components
    {
        sampled_values.clear();
        sampled_mean = 0.0;
    }
    
    // Update empirical proposal distribution
    if ( num_tried > wait_before_learning && num_tried < wait_before_using && num_tried % updates_every == 0)
    {
            
        // store values
        sampled_values.push_back( v );

    }
    else if (num_tried == wait_before_using)
    {
            
        // compute the mean
        double u = sampled_values.size();
        sampled_mean = 0.0;
        size_t updates = sampled_values.size();
        for (size_t i=0; i<updates; ++i)
        {
            // we divide directly to avoid overflow
            sampled_mean += sampled_values[i] / u;
        }
        
        // compute the variance
        u = u - 1.0;
        sampled_var = 0.0;
        for (size_t i=0; i<updates; ++i)
        {
            sampled_var += (sampled_values[i]-sampled_mean)*(sampled_values[i]-sampled_mean) / u;
        }
            
    }

    
    
    
    double ln_Hastings_ratio = 0.0;
    
    if ( stored_index == 0 )
    {
        // draw the new value
        if ( num_tried <= wait_before_using )
        {
            d.redrawValueByIndex( 1 );
            
            // get the base distribution
            TypedDistribution<double> &baseDistribution = d.getBaseDistribution();
            
            // store the proposal ratio
            ln_Hastings_ratio = - baseDistribution.computeLnProbability();
        }
        else
        {
            // set index
            d.setCurrentIndex( 1 );
            
            double var = sampled_var;

            // draw value
            if ( proposal_distribution == NORMAL )
            {
                double sd = sqrt(var);
                double rv = RbStatistics::Normal::rv(sampled_mean, sd, *GLOBAL_RNG);
                d.setValue( new double(rv) );
                
                // store the proposal ratio
                ln_Hastings_ratio -= RbStatistics::Normal::lnPdf(sampled_mean, sd, rv);
            }
            else
            {
                throw RbException("You have selected a proposal distribution for the adaptive-RJ-switch move that is currently not implemented.");
            }
        }
    }
    else
    {
        
        if ( num_tried <= wait_before_using )
        {
            // get the base distribution
            TypedDistribution<double> &baseDistribution = d.getBaseDistribution();
            
            // store the proposal ratio
            ln_Hastings_ratio = baseDistribution.computeLnProbability();
        }
        else
        {
            
            double var = sampled_var;

            // draw value
            if ( proposal_distribution == NORMAL )
            {
                double sd = sqrt(var);
                
                // store the proposal ratio
                ln_Hastings_ratio += RbStatistics::Normal::lnPdf(sampled_mean, sd, v);
            }
        }
        // draw the new value
        d.redrawValueByIndex( 0 );
    }
    
    return ln_Hastings_ratio;
}


/**
 *
 */
void RevBayesCore::AdaptiveReversibleJumpProposal::prepareProposal( void )
{
    
}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void RevBayesCore::AdaptiveReversibleJumpProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "mean";
    if (name_only == false)
    {
        o << " = " << sampled_mean;
    }
    o << ", var";
    if (name_only == false)
    {
        o << " = " << sampled_var;
    }
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void RevBayesCore::AdaptiveReversibleJumpProposal::undoProposal( void )
{
    
    // swap current value and stored value
    variable->setValue( new double(stored_value) );
    
    // also reset the index
    ReversibleJumpMixtureConstantDistribution<double> &d = static_cast< ReversibleJumpMixtureConstantDistribution<double>& >( variable->getDistribution() );
    d.setCurrentIndex( stored_index );
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void RevBayesCore::AdaptiveReversibleJumpProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    variable = static_cast<StochasticNode<double>* >(newN) ;
    
}


void RevBayesCore::AdaptiveReversibleJumpProposal::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void RevBayesCore::AdaptiveReversibleJumpProposal::tune( double rate )
{
    // nothing to do here.
    
}
