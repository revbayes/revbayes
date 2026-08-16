#include "StochasticNodeBase.h"

#include "DagNode.h"
#include "Distribution.h"
#include "DynamicNodeBase.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"

#include <cassert>
#include <cmath>

using namespace RevBayesCore;

/** Take ownership of a distribution and initialize stochastic bookkeeping. */
StochasticNodeBase::StochasticNodeBase(Distribution *d) :
    clamped( false ),
    ignore_data( false ),
    ignore_redraw( false ),
    integrated_out( false ),
    lnProb(),
    stored_ln_prob(),
    distribution( d )
{
}


/** Clone the distribution and flags without copying probability snapshots. */
StochasticNodeBase::StochasticNodeBase(const StochasticNodeBase &n) :
    clamped( n.clamped ),
    ignore_data( n.ignore_data ),
    ignore_redraw( n.ignore_redraw ),
    integrated_out( n.integrated_out ),
    lnProb(),
    stored_ln_prob(),
    distribution( n.distribution->clone() )
{
}


/** Destroy the distribution owned by this implementation base. */
StochasticNodeBase::~StochasticNodeBase(void)
{
    delete distribution;
}


/** Set whether the current value is treated as observed data. */
void StochasticNodeBase::setClamped(bool tf)
{
    clamped = tf;
}


/** Retain the current value but stop treating it as observed data. */
void StochasticNodeBase::unclamp(void)
{
    clamped = false;
}


/** Replace the distribution and copied flags while preserving the existing stored snapshot state. */
void StochasticNodeBase::assign(const StochasticNodeBase &n, DagNode &owner)
{
    // Remove us as the child of the distribution parameters
    detachFromDistributionParameters( owner );

    // Delete the distribution
    delete distribution;

    // Recreate the distribution
    distribution = n.distribution->clone();

    // Get the parameters from the new distribution and add us as child of them in the DAG
    attachToDistributionParameters( owner );

    clamped = n.clamped;
    ignore_data = n.ignore_data;
    ignore_redraw = n.ignore_redraw;
    integrated_out = n.integrated_out;
    lnProb = {};
}


/** Register the stochastic node as a child of every distribution parameter. */
void StochasticNodeBase::attachToDistributionParameters(DagNode &owner) const
{
    // Get the parameters from the distribution and add us as a child of them in the DAG
    const std::vector<const DagNode*> &parents = distribution->getParameters();
    for ( std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->addChild( &owner );

        // Increment the reference count
        // We don't want this parent to get deleted while we are still alive
        (*it)->incrementReferenceCount();
    }
}


/** Draw a new value from the distribution and invalidate the node. */
void StochasticNodeBase::bootstrap(DagNode &owner)
{
    distribution->bootstrap();

    // touch this node for probability recalculation
    owner.touch();
}


/** Sum over integrated parent states recursively using a numerically stable log-sum-exp. */
double StochasticNodeBase::computeRecursiveIntegratedLnProbability(RbOrderedSet<DagNode*> &integratedParents, size_t index)
{
    if ( integratedParents.size() <= index )
    {
        return distribution->computeLnProbability();
    }

    DagNode *parent = integratedParents[index];
    size_t numElements = parent->getNumberOfMixtureElements();
    std::vector<double> lnProbs( numElements, 0.0 );
    double maxLnProb = RbConstants::Double::neginf;

    for ( size_t i = 0; i < numElements; ++i )
    {
        parent->setIntegrationIndex( i );
        lnProbs[i] = computeRecursiveIntegratedLnProbability( integratedParents, index + 1 );
        if ( lnProbs[i] > maxLnProb )
        {
            maxLnProb = lnProbs[i];
        }
    }

    std::vector<double> mixtureProbs = parent->getMixtureProbabilities();
    double probability = 0.0;
    for ( size_t i = 0; i < numElements; ++i )
    {
        probability += std::exp( lnProbs[i] - maxLnProb ) * mixtureProbs[i];
    }

    return std::log( probability ) + maxLnProb;
}


/** Remove the stochastic node's incoming edges and release its references to its parameters. */
void StochasticNodeBase::detachFromDistributionParameters(DagNode &owner) const
{
    // Remove us as the child of the distribution parameters
    // Copy the vector because releasing a parameter may delete it and invalidate related storage.
    std::vector<const DagNode*> parents = distribution->getParameters();
    for ( std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->removeChild( &owner );

        // Decrement the reference count and check whether we need to delete the DAG node
        if ( (*it)->decrementReferenceCount() == 0 )
        {
            delete *it;
        }
    }
}


/** Record this likelihood as affected, or pass discovery downstream when integrated out. */
void StochasticNodeBase::getAffected(DagNode &owner, RbOrderedSet<DagNode*> &affected, const DagNode *affecter)
{
    if ( owner.isIntegratedOut() == false )
    {
        // Insert this node as one of the affected
        affected.insert( &owner );

        // Call the distribution for potential specialized handling (e.g. internal flags)
        distribution->getAffected( affected, affecter );
    }
    else
    {
        // Dispatch the touch message to downstream nodes
        owner.getAffectedNodes( affected );
    }
}


/** Collect integrated stochastic parents, recursively including their integrated ancestors. */
void StochasticNodeBase::getIntegratedParents(const DagNode &owner, RbOrderedSet<DagNode*> &integratedParents) const
{
    std::vector<const DagNode*> parents = owner.getParents();

    // delegate up the DAG
    for ( std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        if ( (*it)->isIntegratedOut() == true )
        {
            (*it)->getIntegratedParents( integratedParents );
            integratedParents.insert( const_cast<DagNode*>( *it ) );
        }
    }
}


/** Compute and cache this node's current log probability when needed. */
double StochasticNodeBase::getLnProbability(DagNode &owner)
{
    if ( not lnProb )
    {
        // compute and store log-probability
        if ( integrated_out || ignore_data )
        {
            lnProb = 0.0;
        }
        else
        {
            RbOrderedSet<DagNode*> integratedParents;
            owner.getIntegratedParents( integratedParents );
            lnProb = computeRecursiveIntegratedLnProbability( integratedParents, 0 );
        }
    }

    return lnProb.value();
}


/** Compute likelihoods for every mixture state, optionally returning them on the log scale. */
std::vector<double> StochasticNodeBase::getMixtureLikelihoods(const DagNode &owner, bool useLog) const
{
    std::vector<double> lnProbs;
    if ( owner.isIntegratedOut() == false )
    {
        // TODO: This is only for safety. We should actually never get in here!
        lnProbs.push_back( distribution->computeLnProbability() );
        return lnProbs;
    }

    // temporarily disable integration
    integrated_out = false;

    size_t numElements = owner.getNumberOfMixtureElements();
    lnProbs.assign( numElements, 0.0 );
    std::vector<double> probabilities( numElements, 0.0 );
    double maxLnProb = RbConstants::Double::neginf;
    RbOrderedSet<DagNode*> affected;
    owner.getAffectedNodes( affected );

    for ( size_t i = 0; i < numElements; ++i )
    {
        const_cast<DagNode&>( owner ).setIntegrationIndex( i );
        for ( size_t j = 0; j < affected.size(); ++j )
        {
            lnProbs[i] += affected[j]->getLnProbability();
        }
        if ( RbMath::isNan( lnProbs[i] ) == true )
        {
            lnProbs[i] = RbConstants::Double::neginf;
        }
        if ( lnProbs[i] > maxLnProb )
        {
            maxLnProb = lnProbs[i];
        }
    }

    std::vector<double> mixtureProbs = owner.getMixtureProbabilities();
    double probability = 0.0;
    for ( size_t i = 0; i < numElements; ++i )
    {
        probability += std::exp( lnProbs[i] - maxLnProb ) * mixtureProbs[i];
    }

    double lnProbability = std::log( probability ) + maxLnProb;
    for ( size_t i = 0; i < numElements; ++i )
    {
        probabilities[i] = std::exp( lnProbs[i] - maxLnProb ) * mixtureProbs[i] / probability;
        lnProbs[i] += std::log( mixtureProbs[i] ) - lnProbability;
    }

    // switch back integration flag
    integrated_out = true;
    return useLog ? lnProbs : probabilities;
}


/** Return the mixture probabilities supplied by the distribution. */
std::vector<double> StochasticNodeBase::getMixtureProbabilities(void) const
{
    return distribution->getMixtureProbabilities();
}


/** Return the number of mixture values supplied by the distribution. */
size_t StochasticNodeBase::getNumberOfMixtureElements(void) const
{
    return distribution->getNumberOfMixtureElements();
}


/**
 * Get the parents of this node. Simply ask the distribution to provide its parameters,
 * no need to keep parents here.
 */
std::vector<const DagNode*> StochasticNodeBase::getParents(void) const
{
    return distribution->getParameters();
}


/** Return the stored log probability, distinguishing absent and never-computed snapshots. */
double StochasticNodeBase::getPrevLnProbability(void) const
{
    /*
     * NOTE: If there is no previous probability then we could do a few things:
     *         (1) throw an exception (current done).
     *         (2) return an optional<double> to indicate if there is a previous probability or not.
     *         (3) return the current probability. This would make the method non-const.
     *       Right now we do (1).
     */

    // 1. If the node is not affected/touched, then throw an exception.
    if ( not stored_ln_prob )
    {
        throw RbException() << "getPrevLnProbability: no previous probability!";
    }
    // 2. If we touched the node when the log probability was not calculated, then we don't have a value for
    // the probability of the previous state.
    if ( not *stored_ln_prob )
    {
        throw RbException() << "getLnProbabilityRatio: the log probability for the previous state was never calculated";
    }

    // 3. If (a) the node is touched/affected and (b) we know the previous probability, then use it.
    return **stored_ln_prob;
}


/** Report whether the node's current value is observed data. */
bool StochasticNodeBase::isClamped(void) const
{
    return clamped;
}


/** Report whether the node's likelihood is currently ignored. */
bool StochasticNodeBase::isIgnoredData(void) const
{
    return ignore_data;
}


/** Report whether the node is marginalized over its mixture values. */
bool StochasticNodeBase::isIntegratedOut(void) const
{
    return integrated_out;
}


/**
 * Keep the current value of the node.
 * At this point, we also need to make sure we update the stored ln probability.
 */
void StochasticNodeBase::keepMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *affecter)
{
    if ( dynamicState.isTouched() == true )
    {
        stored_ln_prob = {};

        if ( not lnProb )
        {
            if ( integrated_out || ignore_data )
            {
                lnProb = 0.0;
            }
            else
            {
                RbOrderedSet<DagNode*> integratedParents;
                owner.getIntegratedParents( integratedParents );
                lnProb = computeRecursiveIntegratedLnProbability( integratedParents, 0 );
            }
        }

        distribution->keep( affecter );

        // clear the list of touched element indices
        owner.clearTouchedElementIndices();
        if ( owner.isIntegratedOut() == true )
        {
            // Dispatch the touch message to downstream nodes
            owner.keepAffected();
        }
    }

    assert( lnProb );

    // delegate call
    dynamicState.keepMe();
}


/** Redraw an unclamped node unless redraws are suppressed, then invalidate its probability. */
void StochasticNodeBase::redraw(DagNode &owner, SimulationCondition condition)
{
    // draw the value
    if ( ignore_redraw == false )
    {
        if ( owner.isClamped() )
        {
            throw RbException( "Cannot modify the value of a clamped node." );
        }
        distribution->redrawValue( condition );
    }

    // touch this node for probability recalculation
    owner.touch();
}


/** Notify the distribution that its model has been reinitialized. */
void StochasticNodeBase::reInitializeMe(void)
{
    distribution->reInitialized();
}


/** Restore the old value of the node and tell affected. */
void StochasticNodeBase::restoreMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *restorer)
{
    if ( dynamicState.isTouched() == true )
    {
        lnProb = stored_ln_prob.value();
        stored_ln_prob = {};    // An almost impossible value for the density

        // reset flags that recalculation is not needed
        assert( lnProb );

        // call for potential specialized handling (e.g. internal flags)
        distribution->restore( restorer );

        // clear the list of touched element indices
        owner.clearTouchedElementIndices();
        if ( owner.isIntegratedOut() == true )
        {
            // Dispatch the touch message to downstream nodes
            owner.restoreAffected();
        }
    }

    // delegate call
    dynamicState.restoreMe();
}


/** Set the active PID of this specific DAG node object. */
void StochasticNodeBase::setActivePIDSpecialized(size_t activePid, size_t numProcesses)
{
    if ( distribution != NULL )
    {
        distribution->setActivePID( activePid, numProcesses );
    }
}


/** Ignore a clamped node's likelihood and reset its probability cache to zero. */
void StochasticNodeBase::setIgnoreData(const DagNode &owner, bool tf)
{
    if ( tf && owner.isClamped() == false )
    {
        throw RbException() << "Error: cannot ignore data at node '" << owner.getName() << "' because it is not clamped (has no data)!";
    }

//    PROBLEM: Right now vectors are children of their elements, so this prohibits using vectors.
//             We check for stochastic descendants?
//    if (tf and not this->children.empty())
//        throw RbException()<<"Error: cannot ignore data at node '"<<this->getName()<<"' because it has children! (e.g. "<<this->children[0]->getName()<<")";

    ignore_data = tf;
    assert( not stored_ln_prob );
    lnProb = 0.0;
}


/** Set whether redraw requests should retain the current value. */
void StochasticNodeBase::setIgnoreRedraw(bool tf)
{
    ignore_redraw = tf;
}


/** Set whether downstream likelihoods marginalize over this node. */
void StochasticNodeBase::setIntegratedOut(bool tf)
{
    integrated_out = tf;
}


/** Propagate the requested MCMC mode to the distribution. */
void StochasticNodeBase::setMcmcMode(bool tf)
{
    distribution->setMcmcMode( tf );
}


/**
 * This function replaces the earlier swapParameter function. If we rely on the
 * internal RevBayesCore::Distribution to manage our parents, we simply need to ask
 * the distribution to swap its parameter, and then manage the connection of the
 * old and new parents (parameters) to this node.
 */
void StochasticNodeBase::swapParent(DagNode &owner, const DagNode *oldParent, const DagNode *newParent)
{
    // We are sure to get into trouble if either one of these is NULL
    if ( oldParent == NULL || newParent == NULL )
    {
        throw RbException( "Attempt to swap NULL distribution parameter of RevBayesCore::StochasticNode" );
    }

    // This throws an error if the oldParent cannot be found
    distribution->swapParameter( oldParent, newParent );

    oldParent->removeChild( &owner );
    if ( oldParent->decrementReferenceCount() == 0 )
    {
        delete oldParent;
    }

    newParent->addChild( &owner );
    newParent->incrementReferenceCount();
    owner.touch();
}


/** Touch this node for recalculation. */
void StochasticNodeBase::touchMe(DagNode &owner, DynamicNodeBase &dynamicState, const DagNode *toucher, bool touchAll)
{
    if ( dynamicState.isTouched() == false )
    {
        assert( not stored_ln_prob );
        stored_ln_prob = lnProb;
    }

    lnProb = {};

    // call for potential specialized handling (e.g. internal flags), we might have been touched already by someone else, so we need to delegate regardless
    distribution->touch( toucher, touchAll );

    // delegate call
    dynamicState.touchMe();

    if ( owner.isIntegratedOut() == true )
    {
        // Dispatch the touch message to downstream nodes
        owner.touchAffected( touchAll );
    }
}
