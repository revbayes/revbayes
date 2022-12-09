#include "AdaptiveReversibleJumpProposal.h"
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
AdaptiveReversibleJumpProposal::AdaptiveReversibleJumpProposal( StochasticNode<double> *n, size_t n0, size_t c0, size_t m ) : Proposal(),
    variable( n ),
    stored_value( NULL ),
    stored_index( 0 ),
    wait_before_learning( n0 ),
    wait_before_using ( c0 ),
    max_updates ( m ),
    num_tried ( 0 ),
    updates ( 0 )
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
    max_updates ( p.max_updates ),
    num_tried ( p.num_tried ),
    updates ( p.updates )
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
    
    if ( num_tried == 1 )
    // First time using move, setting up components
    {
        sampled_values.clear()
    }
    else
    // Update empirical proposal distribution
    // However, we only change the matrix being used when we tune the variance parameter
    {
        // Store values
        for ( size_t i=0; i<dim; ++i )
        {
            storedValues[i] = x[i];
            storedValuesUntransformed[i] = x_untransformed[i];
        }

        if ( num_tried > wait_before_learning && updates <= max_updates)
        {
            ++updates;

            double n = double(updates);

            for ( size_t i=0; i<dim; ++i )
            {
                // Update covariances first (to save us tracking current and previous averages)
                for (size_t j=i; j < dim; ++j)
                {
                    C_emp[i][j] = 1/n * ( C_emp[i][j] * (n - 1.0) + ((n - 1.0)/n) * ((x[i] - x_bar[i]) * (x[j] - x_bar[j])) );
                    C_emp[j][i] = C_emp[i][j];
                }

                // Update averages
                x_bar[i] = 1/n * x[i] + (n - 1)/n * x_bar[i];
            }

        }

    }

    // Move
    std::vector<double> x_prime = rMVNCholesky(x, AVMVN_cholesky_L, *rng, sigma);

    // This also sets all x to x_prime
    calculateHastingsRatio(x_prime, x);

    return lnHastingsratio;
    
    double lnHastingsratio = 0.0;
    
    if ( stored_index == 0 )
    {
        // draw the new value
        d.redrawValueByIndex( 1 );

        // get the base distribution
        TypedDistribution<double> &baseDistribution = d.getBaseDistribution();
        
        // store the proposal ratio
        lnHastingsratio = - baseDistribution.computeLnProbability();
    }
    else
    {
        
        // get the base distribution
        TypedDistribution<double> &baseDistribution = d.getBaseDistribution();
    
        // store the proposal ratio
        lnHastingsratio = baseDistribution.computeLnProbability();
        
        // draw the new value
        d.redrawValueByIndex( 0 );
    }
    
    return lnHastingsratio;
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
    // nothing to print
    
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
    variable->setValue( stored_value );
    
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
