#ifndef DeterministicNode_H
#define DeterministicNode_H

#include "DeterministicNodeBase.h"
#include "DynamicNode.h"
#include "FunctionTypeUtilities.h"
#include "RbException.h"
#include "TypedFunction.h"

namespace RevBayesCore {

    template<class valueType>
    class DeterministicNode : public DynamicNode<valueType>, protected DeterministicNodeBase {

    public:
        DeterministicNode(const std::string &n, TypedFunction<valueType> *f);
        DeterministicNode(const DeterministicNode<valueType> &n);                                                                       //!< Copy constructor
        virtual                                            ~DeterministicNode(void);                                                    //!< Virtual destructor

        DeterministicNode&                                 operator=(const DeterministicNode &n);                                       //!< Assignment operator

        // public methods
        void                                               bootstrap(void);                                                             //!< Bootstrap the current value of the node (applies only to stochastic nodes)
        virtual DeterministicNode<valueType>*              clone(void) const;
        virtual TypedFunction<valueType>&                  getFunction(void);
        virtual const TypedFunction<valueType>&            getFunction(void) const;
        void                                               getIntegratedParents(RbOrderedSet<DagNode*> &ip) const;
        double                                             getLnProbability(void);
        double                                             getLnProbabilityRatio(void);
        valueType&                                         getValue(void);
        const valueType&                                   getValue(void) const;
        bool                                               isConstant(void) const;                                                      //!< Is this DAG node constant?
        virtual void                                       printStructureInfo(std::ostream &o, bool verbose=false) const;               //!< Print the structural information (e.g. name, value-type, distribution/function, children, parents, etc.)
        void                                               redraw(SimulationCondition c = SimulationCondition::MCMC);
        void                                               reInitializeMe(void);                                                        //!< The DAG was re-initialized so maybe you want to reset some stuff (delegate to distribution)
        void                                               setMcmcMode(bool tf);                                                        //!< Set the modus of the DAG node to MCMC mode.
        void                                               setValueFromFile(const path &dir);                                           //!< Set value from string.
        void                                               setValueFromString(const std::string &v);                                    //!< Set value from string.

        // Parent DAG nodes management functions
        virtual std::vector<const DagNode*>                getParents(void) const;                                                      //!< Get the set of parents
        virtual void                                       swapParent(const DagNode *oldParent, const DagNode *newParent);              //!< Exchange the parent (function parameter)

    protected:
        void                                               getAffected(RbOrderedSet<DagNode*> &affected, const DagNode *affecter);      //!< Mark and get affected nodes
        void                                               keepMe(const DagNode *affecter);                                             //!< Keep value of this and affected nodes
        void                                               restoreMe(const DagNode *restorer);                                          //!< Restore value of this nodes
        void                                               swapParameter(const DagNode *oldParent, const DagNode *newParent);           //!< Swap the parameter of this node (needs overwriting in deterministic and stochastic nodes)
        virtual void                                       touchMe(const DagNode *toucher, bool touchAll);                              //!< Touch myself and tell affected nodes value is reset
    };

}


/** Construct a typed deterministic node and attach it to its function parameters. */
template<class valueType>
RevBayesCore::DeterministicNode<valueType>::DeterministicNode(const std::string &n, TypedFunction<valueType> *f) :
    DynamicNode<valueType>( n ),
    DeterministicNodeBase( f )
{
    this->type = DagNode::DETERMINISTIC;
    this->attachToFunctionParameters( *this );

    // Set us as the DAG node of the function
    f->setDeterministicNode( this );
}


/** Copy a deterministic node, attach its cloned function, and restore the typed back-pointer. */
template<class valueType>
RevBayesCore::DeterministicNode<valueType>::DeterministicNode(const DeterministicNode<valueType> &n) :
    DynamicNode<valueType>( n ),
    DeterministicNodeBase( n )
{
    this->type = DagNode::DETERMINISTIC;
    this->attachToFunctionParameters( *this );

    // Set us as the DAG node of the function
    assumeFunctionReturns<valueType>( this->function )->setDeterministicNode( this );
}


/** Detach the node before the implementation base destroys its function. */
template<class valueType>
RevBayesCore::DeterministicNode<valueType>::~DeterministicNode(void)
{
    this->detachFromFunctionParameters( *this );
}


/** Assignment operator. Make sure we deal with parent nodes correctly here. */
template<class valueType>
RevBayesCore::DeterministicNode<valueType>& RevBayesCore::DeterministicNode<valueType>::operator=(const DeterministicNode<valueType> &n)
{
    if ( &n != this )
    {
        // Call base class assignment operator
        DynamicNode<valueType>::operator=( n );
        DeterministicNodeBase::assign( n, *this );

        // Set us as the DAG node of the new function
        assumeFunctionReturns<valueType>( this->function )->setDeterministicNode( this );
    }

    return *this;
}


/** Deterministic nodes have no bootstrap action. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::bootstrap(void)
{
    // nothing to do
}


/** Clone this typed deterministic node. */
template<class valueType>
RevBayesCore::DeterministicNode<valueType>* RevBayesCore::DeterministicNode<valueType>::clone(void) const
{
    return new DeterministicNode<valueType>( *this );
}


/**
 * Get the affected nodes.
 * This call is started by the parent. We need to delegate this call to all our children.
 */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::getAffected(RbOrderedSet<DagNode*> &affected, const DagNode* /*affecter*/)
{
    DeterministicNodeBase::getAffected( *this, affected );
}


/** Return the owned function through its typed interface. */
template<class valueType>
RevBayesCore::TypedFunction<valueType>& RevBayesCore::DeterministicNode<valueType>::getFunction(void)
{
    return *assumeFunctionReturns<valueType>( this->function );
}


/** Return the owned function through its const typed interface. */
template<class valueType>
const RevBayesCore::TypedFunction<valueType>& RevBayesCore::DeterministicNode<valueType>::getFunction(void) const
{
    return *assumeFunctionReturns<valueType>( this->function );
}


/** Forward integrated-parent discovery to the non-template implementation. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::getIntegratedParents(RbOrderedSet<DagNode*> &integratedParents) const
{
    DeterministicNodeBase::getIntegratedParents( *this, integratedParents );
}


/** Deterministic nodes contribute no direct log probability. */
template<class valueType>
double RevBayesCore::DeterministicNode<valueType>::getLnProbability(void)
{
    return 0.0;
}


/** Deterministic nodes contribute no direct log-probability ratio. */
template<class valueType>
double RevBayesCore::DeterministicNode<valueType>::getLnProbabilityRatio(void)
{
    return 0.0;
}


/** Return the function parameters as the node's parents. */
template<class valueType>
std::vector<const RevBayesCore::DagNode*> RevBayesCore::DeterministicNode<valueType>::getParents(void) const
{
    return DeterministicNodeBase::getParents();
}


/** Lazily update the function and return its typed value. */
template<class valueType>
valueType& RevBayesCore::DeterministicNode<valueType>::getValue(void)
{
    TypedFunction<valueType> *typedFunction = assumeFunctionReturns<valueType>( this->function );

    // lazy evaluation
    if ( this->needs_update == true || this->force_update == true )
    {
        typedFunction->update();
        this->needs_update = false;
    }

    return typedFunction->getValue();
}


/** Lazily update the function and return its typed value through a const node. */
template<class valueType>
const valueType& RevBayesCore::DeterministicNode<valueType>::getValue(void) const
{
    TypedFunction<valueType> *typedFunction = assumeFunctionReturns<valueType>( this->function );

    // lazy evaluation
    if ( this->needs_update == true || this->force_update == true )
    {
        typedFunction->update();
        this->needs_update = false;
    }

    return typedFunction->getValue();
}


/** Report whether every function parameter is constant. */
template<class valueType>
bool RevBayesCore::DeterministicNode<valueType>::isConstant(void) const
{
    return DeterministicNodeBase::isConstant();
}


/** Print struct for user. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::printStructureInfo(std::ostream &o, bool verbose) const
{
    o << "_dagType      = Deterministic node (function)" << std::endl;
    o << "_function     = <" << this->function << ">" << std::endl;
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


/** Deterministic nodes have no redraw action. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::redraw(SimulationCondition /*c*/)
{
    // nothing to do
    // the touch should have called our update
}


/** Forward model reinitialization to the owned function. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::reInitializeMe(void)
{
    DeterministicNodeBase::reInitializeMe();
}


/** Commit deterministic and dynamic state through the implementation bases. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::keepMe(const DagNode *affecter)
{
    DeterministicNodeBase::keepMe( *this, *this, affecter );
}


/** Restore deterministic and dynamic state through the implementation bases. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::restoreMe(const DagNode *restorer)
{
    DeterministicNodeBase::restoreMe( *this, *this, restorer );
}


/** Deterministic nodes do not propagate an MCMC mode to functions. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::setMcmcMode(bool /*tf*/)
{
    // nothing to do
}


/** Reject attempts to deserialize a deterministic value from a file. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::setValueFromFile(const RevBayesCore::path & /*dir*/)
{
    throw RbException( "Cannot set a deterministic node from a file." );
}


/** Reject attempts to deserialize a deterministic value from a string. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::setValueFromString(const std::string & /*v*/)
{
    throw RbException( "Cannot set a deterministic node from a string." );
}


/** Forward function-parameter replacement to the implementation base. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::swapParent(const DagNode *oldParent, const DagNode *newParent)
{
    DeterministicNodeBase::swapParent( *this, oldParent, newParent );
}


/** Forward deterministic invalidation to the implementation bases. */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::touchMe(const DagNode *toucher, bool touchAll)
{
    DeterministicNodeBase::touchMe( *this, *this, toucher, touchAll );
}

#endif
