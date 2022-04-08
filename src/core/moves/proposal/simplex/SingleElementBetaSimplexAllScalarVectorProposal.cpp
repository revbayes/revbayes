#include "SingleElementBetaSimplexAllScalarVectorProposal.h"

#include <stddef.h>
#include <cmath>

#include "DistributionBeta.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
SingleElementBetaSimplexAllScalarVectorProposal::SingleElementBetaSimplexAllScalarVectorProposal( StochasticNode<Simplex> *n, std::vector<StochasticNode<double> *> sv, double a, double p) : Proposal(p),
simplex ( n ),
scalar_vector( sv ),
alpha( a )
{
    // tell the base class to add the node
    addNode( simplex );
    
    for (std::vector< StochasticNode<double> *>::const_iterator it = scalar_vector.begin(); it != scalar_vector.end(); it++)
    {
        addNode( *it );
    }
    
    stored_scalar_vector = std::vector<double>(scalar_vector.size(), 0.0);
    
}


void SingleElementBetaSimplexAllScalarVectorProposal::addIndex( size_t v )
{
    if (v < simplex->getValue().size())
    {
        indices.insert(v);
    }
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void SingleElementBetaSimplexAllScalarVectorProposal::cleanProposal( void )
{
    if ( failed == false )
    {
        RbOrderedSet<DagNode*> affected;
        simplex->initiatefindUniqueDescendants( affected );
        
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == stored_simplex.size())
            {
                (*it)->clearTouchedElementIndices();
            }
        }
    }
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
SingleElementBetaSimplexAllScalarVectorProposal* SingleElementBetaSimplexAllScalarVectorProposal::clone( void ) const
{
    
    return new SingleElementBetaSimplexAllScalarVectorProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& SingleElementBetaSimplexAllScalarVectorProposal::getProposalName( void ) const
{
    static std::string name = "SingleElementBetaSimplexAllScalarVector";
    
    return name;
}


double SingleElementBetaSimplexAllScalarVectorProposal::getProposalTuningParameter( void ) const
{
    return alpha;
}


/**
 * Perform the proposal.
 *
 * A Beta-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Beta(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
double SingleElementBetaSimplexAllScalarVectorProposal::doProposal( void )
{
    // reset flag
    failed = false;
    
    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;
    
    // now we store all necessary values
    stored_simplex = simplex->getValue();
    const RbVector<double>& simplex_current = simplex->getValue();
    RbVector<double> simplex_new = simplex_current;
    
    size_t chosen_index;
    if ( indices.size() > 0 )
    {
        std::set<size_t>::iterator it = indices.begin();
        if (indices.size() > 1)
        {
            size_t index = size_t( floor(rng->uniform01() * double(indices.size())) );
            std::advance(it, index);
        }
        chosen_index = *it;
    }
    else
    {
        chosen_index = size_t( floor(rng->uniform01() * double(simplex_current.size())) );
    }
    
    double value_current = simplex_current[chosen_index];
    
    // draw new value
    double a = alpha * value_current + 1.0;
    double b = alpha * (1.0 - value_current) + 1.0;
    
    double value_new = RbStatistics::Beta::rv(a, b, *rng);
    simplex_new[chosen_index] = value_new;
    
    double new_a = alpha * value_new + 1.0;
    double new_b = alpha * (1.0 - value_new) + 1.0;
    double forward = RbStatistics::Beta::lnPdf(a, b, value_new);
    double backward = RbStatistics::Beta::lnPdf(new_a, new_b, value_current);
    
    double scaling_factor_unaffected = (1.0 - value_new) / (1.0 - value_current);
    for (size_t i = 0; i < simplex_new.size(); ++i)
    {
        if (i != chosen_index)
        {
            simplex_new[i] *= scaling_factor_unaffected;
        }
        
        if ( simplex_new[i] < 1E-100)
        {
            failed = true;
            return RbConstants::Double::neginf;
        }
    }
    
    simplex->setValue( new Simplex(simplex_new), false);
    
    // compute the Hastings ratio
    double lnHastingsratio = (backward - forward) + log( scaling_factor_unaffected ) * (int(simplex_new.size()) - 2);
    
    for (size_t i = 0; i < scalar_vector.size(); ++i)
    {
        stored_scalar_vector[i] = scalar_vector[i]->getValue();
        double scalar_new = stored_scalar_vector[i] / scaling_factor_unaffected;
        scalar_vector[i]->setValue(new double( scalar_new ), false);
    }
    lnHastingsratio -= log( scaling_factor_unaffected ) * scalar_vector.size();
    
    RbOrderedSet<DagNode*> affected;
    simplex->initiatefindUniqueDescendants( affected );
    
    for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
    {
        if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == simplex_new.size())
        {
            (*it)->addTouchedElementIndex(chosen_index);
        }
    }
    
    return lnHastingsratio;
    
}


/**
 *
 */
void SingleElementBetaSimplexAllScalarVectorProposal::prepareProposal( void )
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
void SingleElementBetaSimplexAllScalarVectorProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
    
    o << "alpha = ";
    if (name_only == false)
    {
        o << alpha;
    }
    
}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void SingleElementBetaSimplexAllScalarVectorProposal::undoProposal( void )
{
    
    // we undo the proposal only if it didn't fail
    if ( failed == false )
    {
        // undo the proposal
        simplex->setValue( new Simplex(stored_simplex), false );
        for (size_t i = 0; i < stored_scalar_vector.size(); ++i)
        {
            scalar_vector[i]->setValue( new double( stored_scalar_vector[i] ), false );
        }
        
        RbOrderedSet<DagNode*> affected;
        simplex->initiatefindUniqueDescendants( affected );
        
        for (RbOrderedSet<DagNode*>::iterator it = affected.begin(); it != affected.end(); ++it)
        {
            if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC && (*it)->getNumberOfElements() == stored_simplex.size())
            {
                (*it)->clearTouchedElementIndices();
            }
        }
    }
    
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void SingleElementBetaSimplexAllScalarVectorProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    
    if (oldN == simplex)
    {
        simplex = static_cast<StochasticNode<Simplex>* >(newN) ;
    }
    else
    {
        for (size_t i = 0; i < scalar_vector.size(); ++i)
        {
            if ( scalar_vector[i] == oldN )
            {
                scalar_vector[i] = static_cast<StochasticNode<double> *>(newN);
            }
        }
    }
    
}


void SingleElementBetaSimplexAllScalarVectorProposal::setProposalTuningParameter(double tp)
{
    alpha = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void SingleElementBetaSimplexAllScalarVectorProposal::tune( double rate )
{
    
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        alpha /= (1.0 + ((rate - p)/(1.0 - p)) );
    }
    else
    {
        alpha *= (2.0 - rate/p);
    }
    
}


