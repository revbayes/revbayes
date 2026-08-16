#ifndef StochasticNode_H
#define StochasticNode_H

#include "Cloner.h"
#include "DistributionTypeUtilities.h"
#include "DynamicNode.h"
#include "IsDerivedFrom.h"
#include "RbException.h"
#include "RbVector.h"
#include "Serializer.h"
#include "StochasticNodeBase.h"
#include "TypedDistribution.h"

#include <cassert>

namespace RevBayesCore {

    template <class valueType>
    class StochasticNode : public DynamicNode<valueType>, public MemberObject<RbVector<double>>, public StochasticNodeBase {

    public:
        StochasticNode(const std::string &n, TypedDistribution<valueType> *d);
        StochasticNode(const StochasticNode<valueType> &n);                                                                             //!< Copy constructor
        virtual                                            ~StochasticNode(void);                                                       //!< Virtual destructor

        // Assignment operator
        StochasticNode&                                    operator=(const StochasticNode &n);                                          //!< Assignment operator

        // Basic utility function
        virtual StochasticNode<valueType>*                 clone(void) const;

        // methods
        void                                               bootstrap(void);                                                             //!< Bootstrap the current value of the node (applies only to stochastic nodes)
        void                                               clamp(valueType *val);                                                       //!< Clamp an observation to this random variable
        void                                               executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const; //!< Map the member methods to internal function calls
        virtual TypedDistribution<valueType>&              getDistribution(void);
        virtual const TypedDistribution<valueType>&        getDistribution(void) const;
        void                                               getIntegratedParents(RbOrderedSet<DagNode*> &ip) const;
        virtual double                                     getLnProbability(void);
        virtual double                                     getPrevLnProbability(void) const;
        virtual double                                     getLnProbabilityRatio(void);
        virtual std::vector<double>                        getMixtureLikelihoods(bool log=true) const;
        virtual std::vector<double>                        getMixtureProbabilities(void) const;
        virtual size_t                                     getNumberOfMixtureElements(void) const;                                     //!< Get the number of elements for this value
        valueType&                                         getValue(void);
        const valueType&                                   getValue(void) const;
        bool                                               isClamped(void) const;                                                       //!< Is this DAG node clamped?
        bool                                               isIntegratedOut(void) const;
        bool                                               isIgnoredData(void) const;
        bool                                               isStochastic(void) const;                                                    //!< Is this DAG node stochastic?
        virtual void                                       printStructureInfo(std::ostream &o, bool verbose=false) const;               //!< Print the structural information (e.g. name, value-type, distribution/function, children, parents, etc.)
        void                                               redraw(SimulationCondition c = SimulationCondition::MCMC);                   //!< Redraw the current value of the node (applies only to stochastic nodes)
        virtual void                                       reInitializeMe(void);                                                        //!< The DAG was re-initialized so maybe you want to reset some stuff (delegate to distribution)
        void                                               setIgnoreRedraw(bool tf=true);
        void                                               setIntegratedOut(bool tf=true);
        virtual void                                       setIntegrationIndex(size_t i);
        void                                               setMcmcMode(bool tf);                                                        //!< Set the modus of the DAG node to MCMC mode.
        virtual void                                       setIgnoreData(bool tf);                                                      //!< Set whether we want to have the probability of the prior only.
        virtual void                                       setValue(valueType *val, bool touch=true);                                   //!< Set the value of this node
        void                                               setValueFromFile(const path &dir);                                           //!< Set value from string.
        void                                               setValueFromString(const std::string &v);                                    //!< Set value from string.

        // Parent DAG nodes management functions
        std::vector<const DagNode*>                        getParents(void) const;                                                      //!< Get the set of parents
        void                                               swapParent(const DagNode *oldParent, const DagNode *newParent);              //!< Exchange the parent (distribution parameter)

    protected:
        virtual double                                     computeRecursiveIntegratedLnProbability(RbOrderedSet<DagNode*> &integratedParents, size_t index);
        virtual void                                       getAffected(RbOrderedSet<DagNode*> &affected, const DagNode *affecter);      //!< Mark and get affected nodes
        virtual void                                       keepMe(const DagNode *affecter);                                             //!< Keep value of this and affected nodes
        virtual void                                       restoreMe(const DagNode *restorer);                                          //!< Restore value of this nodes
        virtual void                                       setActivePIDSpecialized(size_t activePid, size_t numProcesses);              //!< Set the number of processes for this class.
        virtual void                                       touchMe(const DagNode *toucher, bool touchAll);                              //!< Tell affected nodes value is reset
    };

}


/** Construct a typed stochastic node and attach it to its distribution parameters. */
template<class valueType>
RevBayesCore::StochasticNode<valueType>::StochasticNode(const std::string &n, TypedDistribution<valueType> *d)
    : DynamicNode<valueType>( n ),
      StochasticNodeBase( d )
{
    this->type = DagNode::STOCHASTIC;
    this->attachToDistributionParameters( *this );

    // Set us as the DAG node of the distribution
    d->setStochasticNode( this );
}


/** Copy a stochastic node, attach its cloned distribution, and restore the typed back-pointer. */
template<class valueType>
RevBayesCore::StochasticNode<valueType>::StochasticNode(const StochasticNode<valueType> &n)
    : DynamicNode<valueType>( n ),
      StochasticNodeBase( n )
{
    this->type = DagNode::STOCHASTIC;
    this->attachToDistributionParameters( *this );

    // Set us as the DAG node of the distribution
    assumeDistributionOf<valueType>( this->distribution )->setStochasticNode( this );
}


/** Detach the node before the implementation base destroys its distribution. */
template<class valueType>
RevBayesCore::StochasticNode<valueType>::~StochasticNode(void)
{
    this->detachFromDistributionParameters( *this );
}


/** Assignment operator. Make sure we deal with parent nodes correctly here. */
template<class valueType>
RevBayesCore::StochasticNode<valueType>& RevBayesCore::StochasticNode<valueType>::operator=(const StochasticNode<valueType> &n)
{
    if ( &n != this )
    {
        // Call base class assignment operators
        DynamicNode<valueType>::operator=( n );
        StochasticNodeBase::assign( n, *this );

        // Set us as the DAG node of the new distribution
        assumeDistributionOf<valueType>( this->distribution )->setStochasticNode( this );
    }

    return *this;
}


/** Forward bootstrap redraw and invalidation to the implementation base. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::bootstrap(void)
{
    StochasticNodeBase::bootstrap( *this );
}


/** Install an observed typed value and mark the node clamped. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::clamp(valueType *val)
{
    // clamp the node with the value
    // we call set value because some derived classes might have special implementations for setting values (e.g. mixtures)
    setValue( val );
    this->clamped = true;
}


/** Clone this typed stochastic node. */
template<class valueType>
RevBayesCore::StochasticNode<valueType>* RevBayesCore::StochasticNode<valueType>::clone(void) const
{
    return new StochasticNode<valueType>( *this );
}


/** Preserve the virtual recursive integration hook while delegating its default behavior. */
template<class valueType>
double RevBayesCore::StochasticNode<valueType>::computeRecursiveIntegratedLnProbability(RbOrderedSet<DagNode*> &integratedParents, size_t index)
{
    return StochasticNodeBase::computeRecursiveIntegratedLnProbability( integratedParents, index );
}


/** Dispatch the two stochastic-node mixture member methods. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::executeMethod(const std::string &n, const std::vector<const DagNode*> & /*args*/, RbVector<double> &rv) const
{
    if ( n == "lnMixtureLikelihoods" )
    {
        rv = this->getMixtureLikelihoods( true );
    }
    else if ( n == "MixtureLikelihoods" )
    {
        rv = this->getMixtureLikelihoods( false );
    }
    else
    {
        throw RbException() << "A DAG node does not have a member method called '" << n << "'.";
    }
}


/**
 * Get the affected stochastic nodes. We keep track of who issued the call (affecter). The
 * implementation simply inserts this node in the set of affected nodes. The call is not
 * passed on to the children, because the likelihood of descendant stochastic nodes is not
 * affected unless the call comes from this node. In the latter case, the call originates
 * in the getAffectedNodes of this node, and passes on directly to the children
 * without inserting this node in the set of affected nodes. See the DagNode base class for
 * the implementation of getAffectedNodes(...).
 */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::getAffected(RbOrderedSet<DagNode*> &affected, const DagNode *affecter)
{
    StochasticNodeBase::getAffected( *this, affected, affecter );
}


/** Return the owned distribution through its typed interface. */
template<class valueType>
RevBayesCore::TypedDistribution<valueType>& RevBayesCore::StochasticNode<valueType>::getDistribution(void)
{
    return *assumeDistributionOf<valueType>( this->distribution );
}


/** Return the owned distribution through its const typed interface. */
template<class valueType>
const RevBayesCore::TypedDistribution<valueType>& RevBayesCore::StochasticNode<valueType>::getDistribution(void) const
{
    return *static_cast<const TypedDistribution<valueType>*>( this->distribution );
}


/** Forward integrated-parent discovery to the non-template implementation. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::getIntegratedParents(RbOrderedSet<DagNode*> &integratedParents) const
{
    StochasticNodeBase::getIntegratedParents( *this, integratedParents );
}


/** Forward per-mixture likelihood calculation to the non-template implementation. */
template<class valueType>
std::vector<double> RevBayesCore::StochasticNode<valueType>::getMixtureLikelihoods(bool useLog) const
{
    return StochasticNodeBase::getMixtureLikelihoods( *this, useLog );
}


/** Compute or return the cached log probability. */
template<class valueType>
double RevBayesCore::StochasticNode<valueType>::getLnProbability(void)
{
    return StochasticNodeBase::getLnProbability( *this );
}


/** Return the current-to-previous log-probability difference for a touched node. */
template<class valueType>
double RevBayesCore::StochasticNode<valueType>::getLnProbabilityRatio(void)
{
    // 1. If the node is not affected/touched, then the probability is the same for the current and previous state.
    if ( not this->stored_ln_prob )
    {
        return 0.0;
    }
    // 2. If we touched the node when the log probability was not calculated, then we don't have a value for
    // the probability of the previous state.
    if ( not *this->stored_ln_prob )
    {
        throw RbException() << "getLnProbabilityRatio: the log probability for the previous state was never calculated";
    }

    // 3. If (a) the node is touched/affected and (b) we know the previous probability, then use it.
    return getLnProbability() - **this->stored_ln_prob;
}


/** Return the previous log probability from the implementation base. */
template<class valueType>
double RevBayesCore::StochasticNode<valueType>::getPrevLnProbability(void) const
{
    return StochasticNodeBase::getPrevLnProbability();
}


/** Return the distribution's mixture probabilities. */
template<class valueType>
std::vector<double> RevBayesCore::StochasticNode<valueType>::getMixtureProbabilities(void) const
{
    return StochasticNodeBase::getMixtureProbabilities();
}


/** Return the distribution's number of mixture elements. */
template<class valueType>
size_t RevBayesCore::StochasticNode<valueType>::getNumberOfMixtureElements(void) const
{
    return StochasticNodeBase::getNumberOfMixtureElements();
}


/** Return the distribution parameters as the node's parents. */
template<class valueType>
std::vector<const RevBayesCore::DagNode*> RevBayesCore::StochasticNode<valueType>::getParents(void) const
{
    return StochasticNodeBase::getParents();
}


/** Return the mutable typed value owned by the distribution. */
template<class valueType>
valueType& RevBayesCore::StochasticNode<valueType>::getValue(void)
{
    return assumeDistributionOf<valueType>( this->distribution )->getValue();
}


/** Return the const typed value owned by the distribution. */
template<class valueType>
const valueType& RevBayesCore::StochasticNode<valueType>::getValue(void) const
{
    return static_cast<const TypedDistribution<valueType>*>( this->distribution )->getValue();
}


/** Report whether this stochastic node is clamped. */
template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isClamped(void) const
{
    return StochasticNodeBase::isClamped();
}


/** Report whether this stochastic node is integrated out. */
template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isIntegratedOut(void) const
{
    return StochasticNodeBase::isIntegratedOut();
}


/** Report whether this stochastic node's data likelihood is ignored. */
template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isIgnoredData(void) const
{
    return StochasticNodeBase::isIgnoredData();
}


/** Identify this node as stochastic. */
template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isStochastic(void) const
{
    return true;
}


/** Commit stochastic and dynamic state through the implementation bases. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::keepMe(const DagNode *affecter)
{
    StochasticNodeBase::keepMe( *this, *this, affecter );
}


/** Print stochastic-node structure while retaining access to DagNode's protected formatting helpers. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::printStructureInfo(std::ostream &o, bool verbose) const
{
    o << "_dagType      = Stochastic node (distribution)" << std::endl;
    o << "_distribution = <" << this->distribution << ">" << std::endl;
    o << "_clamped      = " << ( this->clamped ? "TRUE" : "FALSE" ) << std::endl;
    o << "_lnProb       = " << const_cast<StochasticNode<valueType>*>( this )->getLnProbability() << std::endl;

    if ( verbose == true )
    {
        o << "_stored_ln_prob = ";
        if ( not this->stored_ln_prob )
            o << "EMPTY";
        else if ( not *this->stored_ln_prob )
            o << "UNCOMPUTED";
        else
            o << **this->stored_ln_prob;
        o << std::endl;
    }

    o << "_parents      = ";
    this->printParents( o, 16, 70, verbose );
    o << std::endl;
    o << "_children     = ";
    this->printChildren( o, 16, 70, verbose );
    o << std::endl;

    if ( verbose == true )
    {
        o << "_dagNode      = " << this->name << " <" << this << ">" << std::endl;
        o << "_refCount     = " << this->getReferenceCount() << std::endl;
        o << "_touched      = " << ( this->touched ? "TRUE" : "FALSE" ) << std::endl;
    }
}


/** Forward redraw behavior to the non-template implementation. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::redraw(SimulationCondition condition)
{
    StochasticNodeBase::redraw( *this, condition );
}


/** Forward model reinitialization to the owned distribution. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::reInitializeMe(void)
{
    StochasticNodeBase::reInitializeMe();
}


/** Restore stochastic and dynamic state through the implementation bases. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::restoreMe(const DagNode *restorer)
{
    StochasticNodeBase::restoreMe( *this, *this, restorer );
}


/** Forward process partitioning to the owned distribution. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setActivePIDSpecialized(size_t activePid, size_t numProcesses)
{
    StochasticNodeBase::setActivePIDSpecialized( activePid, numProcesses );
}


/** Set whether redraw requests should retain the current value. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setIgnoreRedraw(bool tf)
{
    StochasticNodeBase::setIgnoreRedraw( tf );
}


/** Set whether downstream likelihoods marginalize over this node. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setIntegratedOut(bool tf)
{
    StochasticNodeBase::setIntegratedOut( tf );
}


/** Clone and install the selected mixture value. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setIntegrationIndex(size_t i)
{
    TypedDistribution<valueType> *typedDistribution = assumeDistributionOf<valueType>( this->distribution );
    valueType *newValue = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is>::createClone( typedDistribution->getParameterValues()[i] );
    this->setValue( newValue );
}


/** Propagate MCMC mode to the owned distribution. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setMcmcMode(bool tf)
{
    StochasticNodeBase::setMcmcMode( tf );
}


/** Set whether this clamped node's likelihood should be ignored. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setIgnoreData(bool tf)
{
    StochasticNodeBase::setIgnoreData( *this, tf );
}


/** Set the value. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setValue(valueType *val, bool forceTouch)
{
    // set the value
    assumeDistributionOf<valueType>( this->distribution )->setValue( val, true );
    if ( forceTouch == true )
    {
        // touch this node for probability recalculation
        this->touch();
    }
}


/** Read a serialized typed value and install it through the normal value path. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setValueFromFile(const RevBayesCore::path &dir)
{
    Serializer<valueType, IsDerivedFrom<valueType, Serializable>::Is>::ressurectFromFile( &getValue(), dir, this->getName() );

    // delegate to the standard function of setting the value
    this->setValue( &this->getValue() );
}


/** Parse a hidden state or a serialized typed value and invalidate the node. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setValueFromString(const std::string &v)
{
    if ( v.size() >= 2 && v.front() == '\'' && v.back() == '\'' )
    {
        // quoted string — this is a hidden state (e.g. allocation index)
        this->distribution->setHiddenStateFromString( v.substr( 1, v.size() - 2 ) );
        this->touch();
    }
    else
    {
        Serializer<valueType, IsDerivedFrom<valueType, Serializable>::Is>::ressurectFromString( &getValue(), v );
        this->setValue( &this->getValue() );
    }
}


/** Forward distribution-parameter replacement to the implementation base. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::swapParent(const DagNode *oldParent, const DagNode *newParent)
{
    StochasticNodeBase::swapParent( *this, oldParent, newParent );
}


/** Forward stochastic invalidation to the implementation bases. */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::touchMe(const DagNode *toucher, bool touchAll)
{
    StochasticNodeBase::touchMe( *this, *this, toucher, touchAll );
}

#endif
