#ifndef StochasticNode_H
#define StochasticNode_H

#include "DynamicNode.h"
#include "RbException.h"

namespace RevBayesCore {
    
    
    template <class valueType>
    class TypedDistribution;
    
    template <class valueType>
    class StochasticNode : public DynamicNode<valueType>, public MemberObject< RbVector<double> > {
        
    public:
        
        StochasticNode(const std::string &n, TypedDistribution<valueType> *d);
        StochasticNode(const StochasticNode<valueType> &n);                                                                             //!< Copy constructor
        
        virtual                                            ~StochasticNode(void);                                                       //!< Virtual destructor
        
        // Assignment operator
        StochasticNode&                                     operator=(const StochasticNode& n);                                         //!< Assignment operator
        
        // Basic utility function
        virtual StochasticNode<valueType>*                  clone(void) const;
        
        // methods
        void                                                bootstrap(void);                                                            //!< Bootstrap the current value of the node (applies only to stochastic nodes)
        void                                                clamp(valueType *val);                                                      //!< Clamp an observation to this random variable
        void                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const; //!< Map the member methods to internal function calls
        virtual TypedDistribution<valueType>&               getDistribution(void);
        virtual const TypedDistribution<valueType>&         getDistribution(void) const;
        void                                                getIntegratedParents(RbOrderedSet<DagNode *>& ip) const;
        virtual double                                      getLnProbability(void);
        virtual double                                      getPrevLnProbability(void) const;
        virtual double                                      getLnProbabilityRatio(void);
        virtual std::vector<double>                         getMixtureLikelihoods(bool log=true) const;
        virtual std::vector<double>                         getMixtureProbabilities(void) const;
        virtual size_t                                      getNumberOfMixtureElements(void) const;                                                        //!< Get the number of elements for this value
        valueType&                                          getValue(void);
        const valueType&                                    getValue(void) const;
        bool                                                isClamped(void) const;                                                      //!< Is this DAG node clamped?
        bool                                                isIntegratedOut(void) const;
        bool                                                isIgnoredData(void) const;
        bool                                                isStochastic(void) const;                                                   //!< Is this DAG node stochastic?
        virtual void                                        printStructureInfo(std::ostream &o, bool verbose=false) const;              //!< Print the structural information (e.g. name, value-type, distribution/function, children, parents, etc.)
        void                                                redraw(SimulationCondition c = SimulationCondition::MCMC);                  //!< Redraw the current value of the node (applies only to stochastic nodes)
        virtual void                                        reInitializeMe(void);                                                       //!< The DAG was re-initialized so maybe you want to reset some stuff (delegate to distribution)
        virtual void                                        setClamped(bool tf);                                                        //!< Set directly the flag whether this node is clamped.
        void                                                setIgnoreRedraw(bool tf=true);
        void                                                setIntegratedOut(bool tf=true);
        virtual void                                        setIntegrationIndex( size_t i );
        void                                                setMcmcMode(bool tf);                                                       //!< Set the modus of the DAG node to MCMC mode.
        virtual void                                        setIgnoreData(bool tf);                                                 //!< Set whether we want to have the probability of the prior only.
        virtual void                                        setValue(valueType *val, bool touch=true);                                  //!< Set the value of this node
        void                                                setValueFromFile(const path &dir);                                          //!< Set value from string.
        void                                                setValueFromString(const std::string &v);                                   //!< Set value from string.
        void                                                unclamp(void);                                                              //!< Unclamp the variable
        
        // Parent DAG nodes management functions
        std::vector<const DagNode*>                         getParents(void) const;                                                     //!< Get the set of parents
        void                                                swapParent(const DagNode *oldP, const DagNode *newP);                       //!< Exchange the parent (distribution parameter)
        
    protected:
        
        virtual double                                      computeRecursiveIntegratedLnProbability(RbOrderedSet<DagNode *>& ig, size_t idx);
        virtual void                                        getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);    //!< Mark and get affected nodes
        virtual void                                        keepMe(const DagNode* affecter);                                            //!< Keep value of this and affected nodes
        virtual void                                        restoreMe(const DagNode *restorer);                                         //!< Restore value of this nodes
        virtual void                                        setActivePIDSpecialized(size_t i, size_t n);                                //!< Set the number of processes for this class.
        virtual void                                        touchMe(const DagNode *toucher, bool touchAll);                             //!< Tell affected nodes value is reset
        
        // protected members
        bool                                                clamped = false;
        bool                                                ignore_data = false;                                                    //!< The PDF for this node is set to 1, removing the effects of the node.  Only for clamped nodes with no children.
        bool                                                ignore_redraw = false;
        mutable bool                                        integrated_out = false;
        std::optional<double>                               lnProb;                                                                     //!< Current log probability, or empty if not computed.
        std::optional<std::optional<double>>                stored_ln_prob;                                                             //!< Previous log probability if (a) there is a previous state and (b) the log probability for it is computed.
        TypedDistribution<valueType>*                       distribution;
        
    };
}


#include "RbConstants.h"
#include "RbOptions.h"
#include "RbMathLogic.h"
#include "TypedDistribution.h"


#include <cmath>


template<class valueType>
RevBayesCore::StochasticNode<valueType>::StochasticNode( const std::string &n, TypedDistribution<valueType> *d )
    : DynamicNode<valueType>( n ),
      distribution( d )
{
    this->type = DagNode::STOCHASTIC;
    
    // Get the parameters from the distribution and add us as a child of them in the DAG
    const std::vector<const DagNode*>& distParents = distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = distParents.begin(); it != distParents.end(); ++it)
    {
        (*it)->addChild( this );
        
        // Increment the reference count
        // We don't want this parent to get deleted while we are still alive
        (*it)->incrementReferenceCount();
    }
    
    // Set us as the DAG node of the distribution
    distribution->setStochasticNode( this );
}



template<class valueType>
RevBayesCore::StochasticNode<valueType>::StochasticNode( const StochasticNode<valueType> &n )
    : DynamicNode<valueType>( n ),
      clamped( n.clamped ),
      ignore_data( n.ignore_data ),
      ignore_redraw( n.ignore_redraw ),
      integrated_out( n.integrated_out ),
      distribution( n.distribution->clone() )
{
    this->type = DagNode::STOCHASTIC;
    
    // Get the parameters from the distribution and add us as a child of them in the DAG
    const std::vector<const DagNode*>& distParents = distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = distParents.begin(); it != distParents.end(); ++it)
    {
        (*it)->addChild( this );
        
        // Increment the reference count
        // We don't want this parent to get deleted while we are still alive
        (*it)->incrementReferenceCount();
    }
    
    // Set us as the DAG node of the distribution
    distribution->setStochasticNode( this );
}


template<class valueType>
RevBayesCore::StochasticNode<valueType>::~StochasticNode( void )
{
    
    // Remove us as the child of the distribution parameters
    std::vector<const DagNode*> distParents = distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = distParents.begin(); it != distParents.end(); ++it)
    {
        (*it)->removeChild( this );
        
        // Decrement the reference count and check whether we need to delete the DAG node
        if ( (*it)->decrementReferenceCount() == 0)
        {
            delete (*it);
        }
        
    }
    
    delete distribution;
    
}


/**
 * Assignment operator. Make sure we deal with parent nodes correctly here.
 */
template<class valueType>
RevBayesCore::StochasticNode<valueType>& RevBayesCore::StochasticNode<valueType>::operator=( const StochasticNode<valueType>& n )
{
    
    if ( &n != this )
    {
        // Call base class assignment operator
        DynamicNode<valueType>::operator=( n );
        
        // Remove us as the child of the distribution parameters
        const std::vector<const DagNode*>& dist_parents = distribution->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = dist_parents.begin(); it != dist_parents.end(); ++it)
        {
            (*it)->removeChild( this );
            
            // Decrement the reference count and check whether we need to delete the DAG node
            if ( (*it)->decrementReferenceCount() == 0)
            {
                delete (*it);
            }
            
        }
        
        // Delete the distribution
        delete distribution;
        
        // Recreate the distribution
        distribution = n.distribution->clone();
        
        // Get the parameters from the new distribution and add us as child of them in the DAG
        const std::vector<const DagNode*>& new_dist_parents = distribution->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = new_dist_parents.begin(); it != new_dist_parents.end(); ++it)
        {
            (*it)->addChild( this );
            
            // Increment the reference count
            // We don't want this parent to get deleted while we are still alive
            (*it)->incrementReferenceCount();
        }
        
        // Set us as the DAG node of the new distribution
        distribution->setStochasticNode( this );
        
        clamped                             = n.clamped;
        ignore_data                     = n.ignore_data;
        ignore_redraw                       = n.ignore_redraw;
        integrated_out                      = n.integrated_out;
        lnProb                              = {};
    }
    
    return *this;
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::bootstrap( void )
{
    
    distribution->bootstrap();
    
    // touch this node for probability recalculation
    this->touch();
    
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::clamp(valueType *val)
{
    // clamp the node with the value
    // we call set value because some derived classes might have special implementations for setting values (e.g. mixtures)
    setValue( val );
    
    clamped = true;
    
}



template<class valueType>
RevBayesCore::StochasticNode<valueType>* RevBayesCore::StochasticNode<valueType>::clone( void ) const
{
    
    return new StochasticNode<valueType>( *this );
}


template<class valueType>
double RevBayesCore::StochasticNode<valueType>::computeRecursiveIntegratedLnProbability(RbOrderedSet<DagNode *>& integrated_parents, size_t index)
{

    double ln_prob = 0;
    
    if ( integrated_parents.size() <= index )
    {
        ln_prob = distribution->computeLnProbability();
    }
    else
    {
        DagNode* this_parent = integrated_parents[index];
        size_t num_mixture_elements = this_parent->getNumberOfMixtureElements();
        std::vector<double> ln_probs = std::vector<double>(num_mixture_elements, 0.0);
        double max_ln_probs = RbConstants::Double::neginf;
        for (size_t i=0; i<num_mixture_elements; ++i)
        {
            this_parent->setIntegrationIndex( i );
            ln_probs[i] = computeRecursiveIntegratedLnProbability( integrated_parents, index+1 );
            if ( ln_probs[i] > max_ln_probs )
            {
                max_ln_probs = ln_probs[i];
            }
        }
        std::vector<double> mixture_probs = this_parent->getMixtureProbabilities();
        double prob = 0.0;
        for (size_t i=0; i<num_mixture_elements; ++i)
        {
            prob += exp(ln_probs[i] - max_ln_probs) * mixture_probs[i];
        }
        ln_prob = log( prob ) + max_ln_probs;
    }
    
    return ln_prob;
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const
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
void RevBayesCore::StochasticNode<valueType>::getAffected( RbOrderedSet<DagNode*>& affected, const DagNode* affecter )
{
    
    if ( isIntegratedOut() == false )
    {
        // Insert this node as one of the affected
        affected.insert( this );
    
        // Call the distribution for potential specialized handling (e.g. internal flags)
        distribution->getAffected( affected, affecter );
    }
    else
    {
        // Dispatch the touch message to downstream nodes
        this->getAffectedNodes( affected );
    }
}


template<class valueType>
RevBayesCore::TypedDistribution<valueType>& RevBayesCore::StochasticNode<valueType>::getDistribution( void )
{
    
    return *distribution;
}


template<class valueType>
const RevBayesCore::TypedDistribution<valueType>& RevBayesCore::StochasticNode<valueType>::getDistribution( void ) const
{
    
    return *distribution;
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::getIntegratedParents(RbOrderedSet<DagNode *>& integrated_parents) const
{
    
    std::vector<const DagNode*> parents = this->getParents();

    // delegate up the DAG
    for (size_t i=0; i<parents.size(); ++i)
    {
        const DagNode *the_parent = parents[i];
        if ( the_parent->isIntegratedOut() == true )
        {
            the_parent->getIntegratedParents( integrated_parents );
            integrated_parents.insert( const_cast< DagNode* >(the_parent) );
        }
        
    }
    
}


template<class valueType>
std::vector<double> RevBayesCore::StochasticNode<valueType>::getMixtureLikelihoods( bool use_log ) const
{
      
    std::vector<double> ln_probs;
    
    if ( isIntegratedOut() == false )
    {
        // TODO: This is only for safety. We should actually never get in here!
        ln_probs.push_back( distribution->computeLnProbability() );
    }
    else
    {
        // temporarily disable integration
        integrated_out = false;
        
        size_t num_mixture_elements = this->getNumberOfMixtureElements();
        ln_probs = std::vector<double>(num_mixture_elements, 0.0);
        std::vector<double> probs = std::vector<double>(num_mixture_elements, 0.0);
        double max_ln_probs = RbConstants::Double::neginf;
        RbOrderedSet<DagNode *> affected;
        this->getAffectedNodes( affected );
        for (size_t i=0; i<num_mixture_elements; ++i)
        {
            const_cast< StochasticNode<valueType>* >(this)->setIntegrationIndex( i );
            for ( size_t j=0; j<affected.size(); ++j )
            {
                ln_probs[i] += affected[j]->getLnProbability();
            }
            if ( RbMath::isNan(ln_probs[i]) == true )
            {
                ln_probs[i] = RbConstants::Double::neginf;
            }
            if ( ln_probs[i] > max_ln_probs )
            {
                max_ln_probs = ln_probs[i];
            }
        }
        std::vector<double> mixture_probs = this->getMixtureProbabilities();
        double prob = 0.0;
        for (size_t i=0; i<num_mixture_elements; ++i)
        {
            prob += exp(ln_probs[i] - max_ln_probs) * mixture_probs[i];
        }
        
        double ln_prob = log( prob ) + max_ln_probs;
        for (size_t i=0; i<num_mixture_elements; ++i)
        {
            probs[i] = exp(ln_probs[i] - max_ln_probs) * mixture_probs[i] / prob;
            ln_probs[i] += log(mixture_probs[i]) - ln_prob;
        }
        
        if ( use_log == false )
        {
            ln_probs = probs;
        }
        
        // switch back integration flag
        integrated_out = true;
    }
    
    return ln_probs;
}



template<class valueType>
double RevBayesCore::StochasticNode<valueType>::getLnProbability( void )
{
    
    if ( not lnProb )
    {
        // compute and store log-probability
        if ( integrated_out or ignore_data )
            lnProb = 0.0;
        else
        {
            RbOrderedSet<DagNode *> integrated_parents;
            getIntegratedParents(integrated_parents);
            lnProb = computeRecursiveIntegratedLnProbability(integrated_parents,0);
        }
    }

    return lnProb.value();
}


template<class valueType>
double RevBayesCore::StochasticNode<valueType>::getLnProbabilityRatio( void )
{
    // 1. If the node is not affected/touched, then the probability is the same for the current and previous state.
    if (not stored_ln_prob)
        return 0;

    // 2. If we touched the node when the log probability was not calculated, then we don't have a value for
    // the probability of the previous state.
    if (not *stored_ln_prob)
        throw RbException()<<"getLnProbabilityRatio: the log probability for the previous state was never calculated";

    // 3. If (a) the node is touched/affected and (b) we know the previous probability, then use it.
    return getLnProbability() - **stored_ln_prob;
}


template<class valueType>
double RevBayesCore::StochasticNode<valueType>::getPrevLnProbability( void ) const
{
    /*
     * NOTE: If there is no previous probability then we could do a few things:
     *         (1) throw an exception (current done).
     *         (2) return an optional<double> to indicate if there is a previous probability or not.
     *         (3) return the current probability.  This would make the method non-const.
     *       Right now we do (1).
     */

    // 1. If the node is not affected/touched, then throw an exception.
    if (not stored_ln_prob)
        throw RbException()<<"getPrevLnProbability: no previous probability!";

    // 2. If we touched the node when the log probability was not calculated, then we don't have a value for
    // the probability of the previous state.
    if (not *stored_ln_prob)
        throw RbException()<<"getLnProbabilityRatio: the log probability for the previous state was never calculated";

    // 3. If (a) the node is touched/affected and (b) we know the previous probability, then use it.
    return **stored_ln_prob;
}


template<class valueType>
std::vector<double> RevBayesCore::StochasticNode<valueType>::getMixtureProbabilities(void) const
{
    return distribution->getMixtureProbabilities();
}


template<class valueType>
size_t RevBayesCore::StochasticNode<valueType>::getNumberOfMixtureElements(void) const
{
    return distribution->getNumberOfMixtureElements();
}


/**
 * Get the parents of this node. Simply ask the distribution to provide its parameters,
 * no need to keep parents here.
 */
template<class valueType>
std::vector<const RevBayesCore::DagNode*> RevBayesCore::StochasticNode<valueType>::getParents( void ) const
{
    return distribution->getParameters();
}


template<class valueType>
valueType& RevBayesCore::StochasticNode<valueType>::getValue( void )
{
    
    return distribution->getValue();
}


template<class valueType>
const valueType& RevBayesCore::StochasticNode<valueType>::getValue( void ) const
{
    
    return distribution->getValue();
}


template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isClamped( void ) const
{
    
    return clamped;
}


template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isIntegratedOut(void) const
{
    
    return integrated_out;
}


template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isIgnoredData(void) const
{

    return ignore_data;
}


template<class valueType>
bool RevBayesCore::StochasticNode<valueType>::isStochastic( void ) const
{
    
    return true;
}


/**
 * Keep the current value of the node.
 * At this point, we also need to make sure we update the stored ln probability.
 */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::keepMe( const DagNode* affecter )
{
    
    if ( this->touched == true )
    {
        stored_ln_prob = {};

        if ( not lnProb )
        {
            if (integrated_out or ignore_data)
                lnProb = 0.0;
            else
            {
                RbOrderedSet<DagNode *> integrated_parents;
                getIntegratedParents(integrated_parents);
                lnProb = computeRecursiveIntegratedLnProbability(integrated_parents,0);
            }
        }
        
        distribution->keep( affecter );
        
        // clear the list of touched element indices
        this->touched_elements.clear();
        
        if ( isIntegratedOut() == true )
        {
            // Dispatch the touch message to downstream nodes
            this->keepAffected();
        }
        
    }

    assert( lnProb );
    
    
    // delegate call
    DynamicNode<valueType>::keepMe( affecter );
    
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::printStructureInfo( std::ostream &o, bool verbose ) const
{
        
    o << "_dagType      = Stochastic node (distribution)" << std::endl;
    o << "_distribution = " << "<" << distribution << ">" << std::endl;
    o << "_clamped      = " << ( this->clamped ? "TRUE" : "FALSE" ) << std::endl;
    o << "_lnProb       = " << const_cast< StochasticNode<valueType>* >( this )->getLnProbability() << std::endl;
    
    if ( verbose == true)
    {
        o << "_stored_ln_prob = ";
        if (not stored_ln_prob)
            o<< "EMPTY";
        else if (not *stored_ln_prob)
            o<< "UNCOMPUTED";
        else
            o<<**stored_ln_prob;

        o<< std::endl;
    }
    o << "_parents      = ";
    this->printParents(o, 16, 70, verbose);
    o << std::endl;
    
    o << "_children     = ";
    this->printChildren(o, 16, 70, verbose);
    o << std::endl;
    
    if ( verbose == true )
    {
        o << "_dagNode      = " << this->name << " <" << this << ">" << std::endl;
        o << "_refCount     = " << this->getReferenceCount() << std::endl;
        o << "_touched      = " << ( this->touched ? "TRUE" : "FALSE" ) << std::endl;
    }
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::redraw( SimulationCondition c )
{
    
    // draw the value
    if ( ignore_redraw == false )
    {
        distribution->redrawValue( c );
    }
    
    // touch this node for probability recalculation
    this->touch();
    
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::reInitializeMe( void )
{
    
    distribution->reInitialized();
    
}


/** Restore the old value of the node and tell affected */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::restoreMe( const DagNode *restorer )
{
    
    if ( this->touched == true )
    {
        lnProb              = stored_ln_prob.value();
        stored_ln_prob      = {};    // An almost impossible value for the density

        // reset flags that recalculation is not needed
        assert(lnProb);

        // call for potential specialized handling (e.g. internal flags)
        distribution->restore(restorer);

        // clear the list of touched element indices
        this->touched_elements.clear();

        if ( isIntegratedOut() == true )
        {
            // Dispatch the touch message to downstream nodes
            this->restoreAffected();
        }
    }
    
    // delegate call
    DynamicNode<valueType>::restoreMe( restorer );
    
}


/**
 * Set the active PID of this specific DAG node object.
 */
template <class valueType>
void RevBayesCore::StochasticNode<valueType>::setActivePIDSpecialized(size_t a, size_t n)
{
    
    if ( distribution != NULL )
    {
        distribution->setActivePID( a, n );
    }
    
}


/**
 * Set directly the flag whether this node is clamped.
 * The caller needs to be responsible enough to know that we will assume
 * that the current value is the observed value.
 * We could use instead as well a call: clamp( getValue() );
 */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setClamped(bool tf)
{
    
    clamped = tf;
    
}


template <class valueType>
void RevBayesCore::StochasticNode<valueType>::setIgnoreRedraw( bool tf )
{
    
    ignore_redraw = tf;

}


template <class valueType>
void RevBayesCore::StochasticNode<valueType>::setIntegratedOut( bool tf )
{
    
    integrated_out = tf;

}


template <class valueType>
void RevBayesCore::StochasticNode<valueType>::setIntegrationIndex( size_t i )
{
    valueType *new_val = Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( distribution->getParameterValues()[i] );
    this->setValue( new_val );

}



template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setMcmcMode(bool tf)
{
    
    distribution->setMcmcMode( tf );
    
}

template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setIgnoreData(bool tf)
{
    if (tf and not isClamped())
        throw RbException()<<"Error: cannot ignore data at node '"<<this->getName()<<"' because it is not clamped (has no data)!";

    if (tf and not this->children.empty())
        throw RbException()<<"Error: cannot ignore data at node '"<<this->getName()<<"' because it has children! (e.g. "<<this->children[0]->getName()<<")";

    ignore_data = tf;
}

/**
 * Set the value.
 */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setValue(valueType *val, bool force_touch)
{
    // set the value
    distribution->setValue( val, true );
    
    if ( force_touch == true )
    {
        // touch this node for probability recalculation
        this->touch();
    }
    
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setValueFromFile(const RevBayesCore::path &dir)
{
    
    Serializer<valueType, IsDerivedFrom<valueType, RevBayesCore::Serializable>::Is >::ressurectFromFile( &getValue(), dir, this->getName() );
    
    // delegate to the standard function of setting the value
    this->setValue( &this->getValue() );
    
}


template<class valueType>
void RevBayesCore::StochasticNode<valueType>::setValueFromString(const std::string &v)
{
    
    Serializer<valueType, IsDerivedFrom<valueType, RevBayesCore::Serializable>::Is >::ressurectFromString( &getValue(), v );
    
    // delegate to the standard function of setting the value
    this->setValue( &this->getValue() );
    
}


/**
 * This function replaces the earlier swapParameter function. If we rely on the
 * internal RevBayesCore::Distribution to manage our parents, we simply need to ask
 * the distribution to swap its parameter, and then manage the connection of the
 * old and new parents (parameters) to this node.
 */
template <class valueType>
void RevBayesCore::StochasticNode<valueType>::swapParent( const RevBayesCore::DagNode *oldParent, const RevBayesCore::DagNode *newParent )
{
    // We are sure to get into trouble if either one of these is NULL
    if ( oldParent == NULL || newParent == NULL )
    {
        throw RbException( "Attempt to swap NULL distribution parameter of RevBayesCore::StochasticNode" );
    }
    
    // This throws an error if the oldParent cannot be found
    distribution->swapParameter( oldParent, newParent );
    
    oldParent->removeChild( this );
    if ( oldParent->decrementReferenceCount() == 0 )
    {
        delete ( oldParent );
    }
    
    newParent->addChild( this );
    newParent->incrementReferenceCount();
    
    this->touch();
}


/** touch this node for recalculation */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::touchMe( const DagNode *toucher, bool touchAll )
{
    
    if ( this->touched == false )
    {
        assert(not stored_ln_prob);
        stored_ln_prob = lnProb;
    }
    
    lnProb = {};
    
    // call for potential specialized handling (e.g. internal flags), we might have been touched already by someone else, so we need to delegate regardless
    distribution->touch( toucher, touchAll );
    
    // delegate call
    DynamicNode<valueType>::touchMe( toucher, touchAll );
    
    if ( isIntegratedOut() == true )
    {
        // Dispatch the touch message to downstream nodes
        this->touchAffected( touchAll );
    }
    
}


/**
 * Unclamp this node. If I understand this correctly, we
 * can just set the clamped flag to false if we keep the
 * current value. We do not need to tell anyone that we
 * have changed value because we really haven't.
 */
template<class valueType>
void RevBayesCore::StochasticNode<valueType>::unclamp( void )
{
    clamped = false;
}


#endif

