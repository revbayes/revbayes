#ifndef DeterministicNode_H
#define DeterministicNode_H

#include "DynamicNode.h"
#include "TypedFunction.h"

namespace RevBayesCore {

    template<class valueType>
    class DeterministicNode : public DynamicNode<valueType> {

    public:
        DeterministicNode(const std::string &n, TypedFunction<valueType> *f);
        DeterministicNode(const DeterministicNode<valueType> &n);                                                                       //!< Copy constructor
        virtual                                            ~DeterministicNode(void);                                                    //!< Virtual destructor

        DeterministicNode&                                  operator=(const DeterministicNode& n);                                      //!< Assignment operator

        // public methods
        void                                                bootstrap(void);                                                            //!< Bootstrap the current value of the node (applies only to stochastic nodes)
        virtual DeterministicNode<valueType>*               clone(void) const;
        virtual TypedFunction<valueType>&                   getFunction(void);
        virtual const TypedFunction<valueType>&             getFunction(void) const;
        void                                                getIntegratedParents(RbOrderedSet<DagNode *>& ip) const;
        double                                              getLnProbability(void);
        double                                              getLnProbabilityRatio(void);
        valueType&                                          getValue(void);
        const valueType&                                    getValue(void) const;
        bool                                                isConstant(void) const;                                                     //!< Is this DAG node constant?
        virtual void                                        printStructureInfo(std::ostream &o, bool verbose=false) const;              //!< Print the structural information (e.g. name, value-type, distribution/function, children, parents, etc.)
        void                                                redraw(SimulationCondition c = SimulationCondition::MCMC);
        void                                                reInitializeMe(void);                                                       //!< The DAG was re-initialized so maybe you want to reset some stuff (delegate to distribution)
        void                                                setMcmcMode(bool tf);                                                       //!< Set the modus of the DAG node to MCMC mode.
        void                                                setValueFromFile(const path &dir);                                          //!< Set value from string.
        void                                                setValueFromString(const std::string &v);                                   //!< Set value from string.

        // Parent DAG nodes management functions
        virtual std::vector<const DagNode*>                 getParents(void) const;                                                     //!< Get the set of parents
        virtual void                                        swapParent(const DagNode *oldParent, const DagNode *newParent);             //!< Exchange the parent (function parameter)

    protected:
        void                                                getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);    //!< Mark and get affected nodes
        void                                                keepMe(const DagNode* affecter);                                            //!< Keep value of this and affected nodes
        void                                                restoreMe(const DagNode *restorer);                                         //!< Restore value of this nodes
        void                                                swapParameter(const DagNode *oldP, const DagNode *newP);                    //!< Swap the parameter of this node (needs overwriting in deterministic and stochastic nodes)
        virtual void                                        touchMe(const DagNode *toucher, bool touchAll);                             //!< Touch myself and tell affected nodes value is reset

    private:

        TypedFunction<valueType>*                           function;
        mutable bool                                        needs_update;
        bool                                                force_update;
    };

}

#include <cassert>
#include "RbOptions.h"


template<class valueType>
RevBayesCore::DeterministicNode<valueType>::DeterministicNode( const std::string &n, TypedFunction<valueType> *f ) :
    DynamicNode<valueType>( n ),
    function( f ),
    needs_update( true ),
    force_update( f->forceUpdates() )
{
    this->type = DagNode::DETERMINISTIC;

    // Get the parameters from the function and add us as a child of them in the DAG
    const std::vector<const DagNode*>& funcParents = function->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = funcParents.begin(); it != funcParents.end(); ++it)
    {
        (*it)->addChild( this );

        // Increment the reference count
        // We don't want this parent to get deleted while we are still alive
        (*it)->incrementReferenceCount();
    }

    // Set us as the DAG node of the function
    function->setDeterministicNode( this );
}


template<class valueType>
RevBayesCore::DeterministicNode<valueType>::DeterministicNode( const DeterministicNode<valueType> &n ) :
    DynamicNode<valueType>( n ),
    function( n.function->clone() ),
    needs_update( true ),
    force_update( n.function->forceUpdates() )
{
    this->type = DagNode::DETERMINISTIC;

    // Get the parameters from the function and add us as a child of them in the DAG
    const std::vector<const DagNode*>& funcParents = function->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = funcParents.begin(); it != funcParents.end(); ++it)
    {
        (*it)->addChild( this );

        // Increment the reference count
        // We don't want this parent to get deleted while we are still alive
        (*it)->incrementReferenceCount();
    }

    // Set us as the DAG node of the function
    function->setDeterministicNode( this );
}


template<class valueType>
RevBayesCore::DeterministicNode<valueType>::~DeterministicNode( void )
{

    // Remove us as the child of the function parameters
    std::vector<const DagNode*> funcParents = function->getParameters();
    for (std::vector<const DagNode*>::iterator it = funcParents.begin(); it != funcParents.end(); ++it)
    {
        (*it)->removeChild( this );

        // Decrement the reference count and check whether we need to delete the DAG node
        if ( (*it)->decrementReferenceCount() == 0)
        {
            delete (*it);
        }

    }

    // free the memory of the function
    delete function;

}


/**
 * Assignment operator. Make sure we deal with parent nodes correctly here.
 */
template<class valueType>
RevBayesCore::DeterministicNode<valueType>& RevBayesCore::DeterministicNode<valueType>::operator=( const DeterministicNode<valueType>& n )
{

    if ( &n != this )
    {
        // Call base class assignment operator
        DynamicNode<valueType>::operator=( n );

        // Remove us as the child of the function parameters
        const std::set<const DagNode*>& funcParents = function->getParameters();
        for (std::set<const DagNode*>::iterator it = funcParents.begin(); it != funcParents.end(); ++it)
        {
            (*it)->removeChild( this );

            // Decrement the reference count and check whether we need to delete the DAG node
            // The function does not do this for us
            if ( (*it)->decrementReferenceCount() == 0)
                delete (*it);
        }

        // Delete the function
        delete function;

        // Recreate the function
        function = n.function->clone();

        // Get the parameters from the new function and add us as child of them in the DAG
        funcParents = function->getParameters();
        for (std::set<const DagNode*>::iterator it = funcParents.begin(); it != funcParents.end(); ++it)
        {
            (*it)->addChild( this );

            // Increment the reference count
            // We don't want this parent to get deleted while we are still alive
            (*it)->incrementReferenceCount();
        }

        // Set us as the DAG node of the new function
        function->setDeterministicNode( this );

        needs_update = true;
        force_update = function->forceUpdates();
    }

    return *this;
}


template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::bootstrap( void )
{
    // nothing to do
}


template<class valueType>
RevBayesCore::DeterministicNode<valueType>* RevBayesCore::DeterministicNode<valueType>::clone( void ) const
{

    return new DeterministicNode<valueType>( *this );
}

/**
 * Get the affected nodes.
 * This call is started by the parent. We need to delegate this call to all our children.
 */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode* affecter)
{

    this->getAffectedNodes( affected );

}


template<class valueType>
RevBayesCore::TypedFunction<valueType>& RevBayesCore::DeterministicNode<valueType>::getFunction( void )
{

    return *function;
}


template<class valueType>
const RevBayesCore::TypedFunction<valueType>& RevBayesCore::DeterministicNode<valueType>::getFunction( void ) const
{

    return *function;
}


template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::getIntegratedParents(RbOrderedSet<DagNode *>& integrated_parents) const
{
    std::vector<const DagNode*> parents = this->getParents();                                                                     //!< Get the set of parents (empty set here)

    // delegate up the DAG
    for (size_t i=0; i<parents.size(); ++i)
    {
        const DagNode *the_parent = parents[i];
        the_parent->getIntegratedParents( integrated_parents );
    }
}

template<class valueType>
double RevBayesCore::DeterministicNode<valueType>::getLnProbability( void )
{

    return 0.0;
}


template<class valueType>
double RevBayesCore::DeterministicNode<valueType>::getLnProbabilityRatio( void )
{

    return 0.0;
}


/**
 * Get the parents of this node. Simply ask the function to provide its parameters,
 * no need to keep parents here.
 */
template<class valueType>
std::vector<const RevBayesCore::DagNode*> RevBayesCore::DeterministicNode<valueType>::getParents( void ) const
{
    return function->getParameters();
}


template<class valueType>
valueType& RevBayesCore::DeterministicNode<valueType>::getValue( void )
{

    // lazy evaluation
    if ( needs_update == true || force_update == true )
    {
        function->update();
        needs_update = false;
    }

    return function->getValue();
}


template<class valueType>
const valueType& RevBayesCore::DeterministicNode<valueType>::getValue( void ) const
{

    // lazy evaluation
    if ( needs_update == true || force_update == true )
    {
        const_cast<TypedFunction<valueType> *>(function)->update();
        needs_update = false;
    }

    return function->getValue();
}


template<class valueType>
bool RevBayesCore::DeterministicNode<valueType>::isConstant( void ) const
{

    // iterate over all parents and only if all parents are constant then this node is constant too
    const std::vector<const DagNode*>& parents = function->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it)
    {
        if ( (*it)->isConstant() == false )
        {
            return false;
        }
    }

    return true;
}


/**
 * Keep the current value of the node.
 * At this point, we just delegate to the children.
 */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::keepMe( const DagNode* affecter )
{

    // delegate call to base class
    // this will unset the touched flag if it was set
    DynamicNode<valueType>::keepMe( affecter );

    // allow specialized recovery in functions
    function->keep( affecter );

    // delegate call
    this->keepAffected();

    // clear the list of touched element indices
    this->touched_elements.clear();

}


/** Print struct for user */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::printStructureInfo( std::ostream& o, bool verbose ) const
{

    o << "_dagType      = Deterministic node (function)" << std::endl;
    o << "_function     = <" << function << ">" << std::endl;

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
void RevBayesCore::DeterministicNode<valueType>::redraw( SimulationCondition c )
{
    // nothing to do
    // the touch should have called our update
}


template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::reInitializeMe( void )
{

    function->reInitialized();

}


/** Restore the old value of the node and tell affected */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::restoreMe( const DagNode *restorer )
{

    // the value has been changed so we need to flag for recomputing the value
    // we need to do that even if the touched flag is unset because it can already have been unset
    // by a reset call from one of our parameter while another of our parameters wasn't unset
    // that means we need to guarantee that either all of our parameters are restore first (which we cannot guarantee currently)
    // or we need to update our value every time one of our parameters is restored.
    needs_update = true;

    // we just mark ourselves as clean, albeit perhaps not being updated
    DynamicNode<valueType>::restoreMe( restorer );

    // call for potential specialized handling (e.g. internal flags)
    function->restore(restorer);

    // clear the list of touched element indices
    this->touched_elements.clear();

    // delegate call
    this->restoreAffected();

}



template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::setMcmcMode(bool tf)
{

    // nothing to do

}


template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::setValueFromFile(const RevBayesCore::path &dir)
{

    throw RbException("Cannot set a deterministic node from a file.");
}


template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::setValueFromString(const std::string &v)
{

    throw RbException("Cannot set a deterministic node from a string.");
}


/**
 * This function replaces the earlier swapParameter function. If we rely on the
 * internal RevBayesCore::Function to manage our parents, we simply need to ask
 * the function to swap its parameters, and then manage the connection of the
 * parents (parameters) to this node.
 */
template <class valueType>
void RevBayesCore::DeterministicNode<valueType>::swapParent( const RevBayesCore::DagNode *oldParent, const RevBayesCore::DagNode *newParent )
{
    // We are sure to get into trouble if either one of these is NULL
    if ( oldParent == NULL || newParent == NULL )
    {
        throw RbException( "Attempt to swap NULL function parameter of RevBayesCore::DeterministicNode" );
    }

    // This throws an error if the oldParent cannot be found
    function->swapParameter( oldParent, newParent );

    oldParent->removeChild( this );
    if ( oldParent->decrementReferenceCount() == 0 )
    {
        delete ( oldParent );
    }

    newParent->addChild( this );
    newParent->incrementReferenceCount();

    this->touch();
}


/**
 * Touch this node for recalculation.
 *
 *
 */
template<class valueType>
void RevBayesCore::DeterministicNode<valueType>::touchMe( const DagNode *toucher, bool touchAll )
{

    // delegate call to base class
    // this will set the touched flag if it wasn't set already
    DynamicNode<valueType>::touchMe( toucher, touchAll );


    // We need to touch the function always because of specialized touch functionality in some functions, like vector functions.
    // In principle, it would sufficient to do the touch once for each toucher, but we do not keep track of the touchers here.
    // call for potential specialized handling (e.g. internal flags), we might have been touched already by someone else, so we need to delegate regardless
    // This is essential for lazy evaluation
    function->touch( toucher );


    // mark for update
    needs_update = true;

    // Dispatch the touch message to downstream nodes
    this->touchAffected( touchAll );

}


#endif
