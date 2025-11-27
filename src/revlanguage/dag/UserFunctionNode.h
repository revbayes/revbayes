#ifndef UserFunctionNode_H
#define UserFunctionNode_H

#include "DynamicNode.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RbFileManager.h"

#include <vector>

namespace RevLanguage {
    
    /**
     * @brief DAG node for user-defined Rev functions
     *
     * This DAG node class is used for deterministic nodes, where the value is
     * controlled by a user-defined function written in Rev. Unlike other
     * deterministic function nodes used in the program, the UserFunctionNode
     * does not use a RevBayesCore::Function to recompute its value. Instead,
     * it processes the Rev code directly, using the parser code constituting
     * the syntax tree.
     */
    template<typename rlType>
    class UserFunctionNode : public RevBayesCore::DynamicNode<typename rlType::valueType> {
        
    public:
        UserFunctionNode(const std::string& n, UserFunction* fxn);                                                  //!< Basic constructor
        UserFunctionNode(const UserFunctionNode<rlType>& n);                                                        //!< Copy constructor
        
        virtual                                 ~UserFunctionNode(void);                                            //!< Virtual destructor
        
        // Assignment operator
        UserFunctionNode&                       operator=(const UserFunctionNode& x);                               //!< Assignment operator
        
        // Public methods
        void                                    bootstrap(void) {}                                                  //!< Bootstrap (or not)
        UserFunctionNode<rlType>*               clone(void) const;                                                  //!< Type-safe clone
        RevBayesCore::DagNode*                  cloneDAG(RevBayesCore::DagNodeMap &nodesMap, std::map<std::string, const RevBayesCore::DagNode* > &names) const;   //!< Clone the entire DAG connected to this node
        void                                    getIntegratedParents(RevBayesCore::RbOrderedSet<RevBayesCore::DagNode *>& ip) const;
        double                                  getLnProbability(void) { return 0.0; }                              //!< Get ln prob
        double                                  getLnProbabilityRatio(void) { return 0.0; }                         //!< Get ln prob ratio
        typename rlType::valueType&             getValue(void);                                                     //!< Get the value
        const typename rlType::valueType&       getValue(void) const;                                               //!< Get the value (const)
        bool                                    isConstant(void) const;                                             //!< Is this DAG node constant?
        virtual void                            printStructureInfo(std::ostream& o, bool verbose=false) const;      //!< Print structure info
        void                                    redraw(RevBayesCore::SimulationCondition c) {}                      //!< Redraw (or not)
        void                                    setMcmcMode(bool tf);                                               //!< Set the modus of the DAG node to MCMC mode.
        void                                    setValueFromFile(const RevBayesCore::path &fn);                     //!< Set value from string.
        void                                    setValueFromString(const std::string &v);                           //!< Set value from string.
        void                                    update(void);                                                       //!< Update current value
        
        // Parent DAG nodes management functions
        std::vector<const RevBayesCore::DagNode*>  getParents(void) const;                                                                     //!< Get the set of parents
        void                                    swapParent(const RevBayesCore::DagNode *oldParent, const RevBayesCore::DagNode *newParent); //!< Exchange the parent (element variable)
        
    protected:
        void                                    getAffected(RevBayesCore::RbOrderedSet<RevBayesCore::DagNode *>& affected, const RevBayesCore::DagNode* affecter);  //!< Mark and get affected nodes
        void                                    keepMe(const RevBayesCore::DagNode* affecter);                                                    //!< Keep value of this and affected nodes
        void                                    restoreMe(const RevBayesCore::DagNode *restorer);                                                 //!< Restore value of this nodes
        void                                    touchMe(const RevBayesCore::DagNode *toucher, bool touchAll);                                                    //!< Touch myself and tell affected nodes value is reset
        
    private:
        UserFunction*                           userFunction;                                                       //!< The user function used to compute the value
        RevPtr<RevVariable>                     returnVariable;                                                     //!< The current value of the function, wrapped as a variable
    };
    
}


using namespace RevLanguage;


/**
 * Basic constructor of a user function node. We need to add us as a child
 * of the arguments of the function, which we obtain by calling its getParameters
 * function.
 */
template<typename rlType>
UserFunctionNode<rlType>::UserFunctionNode( const std::string& n, UserFunction* fxn ) :
RevBayesCore::DynamicNode<typename rlType::valueType>( n ),
userFunction( fxn ),
returnVariable( NULL )
{
    this->type = RevBayesCore::DagNode::DETERMINISTIC;
    
    // Adopt our parents
    const std::vector<const RevBayesCore::DagNode*>& parents = userFunction->getParameters();
    for ( std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->addChild( this );
        (*it)->incrementReferenceCount();
    }
}


/**
 * Copy constructor. We make a clone of the user function
 * here, and adopt its parents.
 */
template<typename rlType>
UserFunctionNode<rlType>::UserFunctionNode( const UserFunctionNode<rlType>& n ) :
RevBayesCore::DynamicNode<typename rlType::valueType>( n ),
userFunction( n.userFunction->clone() ),
returnVariable( NULL )
{
    this->type = RevBayesCore::DagNode::DETERMINISTIC;
    
    // Adopt our parents
    const std::vector<const RevBayesCore::DagNode*>& parents = userFunction->getParameters();
    for ( std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->addChild( this );
        (*it)->incrementReferenceCount();
    }
}


/**
 * Destructor. We need to delete the user function here and there are
 * also some DAG node management tasks to be performed.
 */
template<typename rlType>
UserFunctionNode<rlType>::~UserFunctionNode( void )
{
    // Detach ourselves from our parents
    const std::vector<const RevBayesCore::DagNode*>& parents = userFunction->getParameters();
    for ( std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        (*it)->removeChild( this );
        if ( (*it)->decrementReferenceCount() == 0 )
            delete (*it);
    }
    
    // Delete the user function
    delete userFunction;
}


/**
 * Assignment operator. We make sure we detach ourselves from
 * the old user function parameters and attach ourselves to the
 * parameters of the new user function. There is no need to update
 * the returnVariable, as std::int64_t as we mark ourselves as dirty.
 */
template<typename rlType>
UserFunctionNode<rlType>& UserFunctionNode<rlType>::operator=( const UserFunctionNode<rlType>& x )
{
    if (this != &x)
    {
        // Detach ourselves from our parents
        const std::vector<const RevBayesCore::DagNode*>& oldParents = userFunction->getParameters();
        for ( std::vector<const RevBayesCore::DagNode*>::const_iterator it = oldParents.begin(); it != oldParents.end(); ++it )
        {
            (*it)->removeChild( this );
            if ( (*it)->decrementReferenceCount() == 0 )
                delete (*it);
        }
        
        // Delete the user function
        delete userFunction;
        
        // Call the base class assignment operator
        RevBayesCore::DynamicNode<typename rlType::valueType>::operator=( x );
        
        // Make a clone of the new function
        userFunction = x.userFunction->clone();
        
        // Adopt our new parents
        const std::vector<const RevBayesCore::DagNode*>& newParents = userFunction->getParameters();
        for ( std::vector<const RevBayesCore::DagNode*>::const_iterator it = newParents.begin(); it != newParents.end(); ++it )
        {
            (*it)->addChild( this );
            (*it)->incrementReferenceCount();
        }
        
        // Tell everybody we have been changed
        this->touch();
    }
    
    return *this;
}


/** Type-safe clone function */
template<typename rlType>
UserFunctionNode<rlType>* UserFunctionNode<rlType>::clone( void ) const
{
    return new UserFunctionNode<rlType>( *this );
}


/**
 * Clone the entire graph: clone children, swap parents. Because the argument DAG nodes are
 * doubly managed by us and by the variable wrapper around the argument variable objects, we
 * need a special version of the cloneDAG function. Otherwise we cannot distinguish between the
 * call to swapParent from cloneDAG (where we need to replace the DAG node inside the
 * variable) and the call from the variable wrapper (where we simply should update the
 * connections, and the variable wrapper replaces the DAG node inside it).
 *
 * @todo For this to be perfectly safe, we need all DAG nodes to be connected to their
 *       Rev object variables. The hack below (two variables sharing the same DAG node)
 *       should be safe for the current implementation of model, however, which does
 *       not expose its model variables to the parser.
 */
template<typename rlType>
RevBayesCore::DagNode* UserFunctionNode<rlType>::cloneDAG( RevBayesCore::DagNodeMap& newNodes, std::map<std::string, const RevBayesCore::DagNode* > &names ) const
{
    // Return our clone if we have already been cloned
    if ( newNodes.find( this ) != newNodes.end() )
        return ( newNodes[ this ] );
    
    // We avoid the conflict between calls to swapParent from replacement of the Rev object in an
    // argument variable, and the need to swap parents while cloning the DAG by implementing our
    // own cloneDAG function, where we use the temporary hack of replacing argument variables with
    // variables having new objects with the required DAG nodes.
    UserFunctionNode<rlType>* copy = this->clone();
    
    // Add this node and its copy to the map
    newNodes[ this ] = copy;
    
    // We need to remove the copy as a child of our parents in order to stop recursive calls to
    // cloneDAG on our copy, its copy, etc, when we call cloneDAG on our parents
    const std::vector<const RevBayesCore::DagNode*>& parents = this->getParents();
    for ( std::vector<const RevBayesCore::DagNode*>::const_iterator i = parents.begin(); i != parents.end(); ++i )
    {
        const RevBayesCore::DagNode* theParam = (*i);
        
        theParam->removeChild( copy );
        theParam->decrementReferenceCount();
    }
    
    // Now replace the parents of the copy (which are now the same as our parents) with the parent clones.
    // Here we use a hack where the copy gets arguments that are wrappers around independent Rev objects
    // that have the appropriate DAG node clones inside them. This should work for now.
    std::vector<Argument>& args     = userFunction->getArguments();
    std::vector<Argument>& copyArgs = copy->userFunction->getArguments();
    for ( size_t i = 0; i < args.size(); ++i )
    {
        // Replace the i-th copy argument with an appropriate object
        
        const RevObject &rObject = args[i].getVariable()->getRevObject();
        if ( rObject.isModelObject() )
        {
            const AbstractModelObject &modelObj = static_cast< const AbstractModelObject &>( rObject );
            RevBayesCore::DagNode* the_node = modelObj.getDagNode();
            std::map<std::string, const RevBayesCore::DagNode*> names;
            AbstractModelObject *rObjClone = modelObj.clone();
            rObjClone->setDagNode( the_node->cloneDAG( newNodes, names ) );
            RevPtr<RevVariable> theArgumentVariableClone = new RevVariable( rObjClone, args[i].getVariable()->getName() );
            copyArgs[i] = Argument( theArgumentVariableClone, args[i].getLabel(), args[i].isConstant() );
        
            // Manage the DAG node connections
            RevBayesCore::DagNode* theArgumentNodeClone = theArgumentVariableClone->getRevObject().getDagNode();
            theArgumentNodeClone->addChild( copy );
            theArgumentNodeClone->incrementReferenceCount();
        }
        else
        {
            throw RbException("Missing implementation");
        }

    }
    
    /* Make sure the children clone themselves */
    for ( std::vector<RevBayesCore::DagNode*>::const_iterator it = this->children.begin(); it != this->children.end(); ++it )
        (*it)->cloneDAG( newNodes, names );
    
    return copy;
}


/**
 * Get the affected nodes. This call is started by a parent. We need to delegate this call to all our children.
 */
template<typename rlType>
void UserFunctionNode<rlType>::getAffected( RevBayesCore::RbOrderedSet<RevBayesCore::DagNode *>& affected, const RevBayesCore::DagNode* affecter )
{
    this->getAffectedNodes( affected );
}

template<typename rlType>
void UserFunctionNode<rlType>::getIntegratedParents(RevBayesCore::RbOrderedSet<RevBayesCore::DagNode *>& ip) const
{
    // not applicable
}



/**
 * Get the parents of this node. We simply return the parameters of the user-defined
 * function.
 */
template<typename rlType>
std::vector<const RevBayesCore::DagNode*> UserFunctionNode<rlType>::getParents( void ) const
{
    // Get the parameters from the user-defined function
    return userFunction->getParameters();
}


/** Get the value (non-const version) */
template<typename rlType>
typename rlType::valueType& UserFunctionNode<rlType>::getValue( void )
{
    if ( this->touched )
        update();
    
    return static_cast< RevBayesCore::TypedDagNode<typename rlType::valueType>* >( returnVariable->getRevObject().getDagNode() )->getValue();
}


/**
 * Get the value (const version). Note that we need to type-cast here to get access to the
 * update function in case we have been touched and have not updated before. This is required
 * by the lazy evaluation mechanism but it is a perfectly safe const cast. If we had not used
 * lazy evaluation but immediate update, the const cast had not been needed.
 */
template<typename rlType>
const typename rlType::valueType& UserFunctionNode<rlType>::getValue( void ) const
{
    if ( this->touched )
        const_cast<UserFunctionNode<rlType>*>( this )->update();
    
    return static_cast< RevBayesCore::TypedDagNode<typename rlType::valueType>* >( returnVariable->getRevObject().getDagNode() )->getValue();
}


/**
 * Is this node constant? We ask our parents (the arguments of the function) if they are all constant.
 * If so, we are also constant, otherwise we are not constant.
 */
template<typename rlType>
bool UserFunctionNode<rlType>::isConstant( void ) const
{
    // Iterate over all parents and only if all parents are constant then this node is constant too
    const std::vector<const RevBayesCore::DagNode*>& parents = getParents();
    for ( std::vector<const RevBayesCore::DagNode*>::const_iterator it = parents.begin(); it != parents.end(); ++it )
    {
        if ( !(*it)->isConstant() )
        {
            return false;
        }
    }
    
    return true;
}


/**
 * Keep the current value of the node. We need not and should not change the touched
 * flag here. If we have not been updated, we should just leave the touched flag in
 * the dirty state. We pass on the message unconditionally, to be on the safe side
 * if the DAG is in an inconsistent state.
 */
template<typename rlType>
void UserFunctionNode<rlType>::keepMe( const RevBayesCore::DagNode* affecter )
{
    
    // Pass the call to downstream nodes
    this->keepAffected();
}


/** Print struct for user */
template< typename rlType >
void UserFunctionNode<rlType>::printStructureInfo( std::ostream& o, bool verbose ) const
{
    if ( verbose )
    {
        o << "_dagNode      = " << this->name << " <" << this << ">" << std::endl;
    }
    else
    {
        if ( this->name != "")
            o << "_dagNode      = " << this->name << std::endl;
        else
            o << "_dagNode      = <" << this << ">" << std::endl;
    }
    
    o << "_dagType      = Deterministic Rev user function node" << std::endl;
    
    if ( verbose == true )
    {
        o << "_refCount     = " << this->getReferenceCount() << std::endl;
    }
    
    o << "_function     = " << userFunction->getRevDeclaration() << std::endl;
    
    
    if ( verbose == true )
    {
        o << "_touched      = " << ( this->touched ? "TRUE" : "FALSE" ) << std::endl;
    }
    
    o << "_parents      = ";
    this->printParents( o, 16, 70, verbose );
    o << std::endl;
    
    o << "_children     = ";
    this->printChildren( o, 16, 70, verbose );
    
    o << std::endl;
}


/**
 * We get this call during the restore phase of an MCMC simulation. We do
 * not know for sure whether the return variable has been updated during
 * the touched phase of the simulation. If it has, it is important to mark
 * us as dirty, but we do not need to update the value immediately. Instead,
 * we rely on lazy evaluation, and update the value when somebody asks for
 * the value.
 *
 * The call is passed on unconditionally, to guard us against inconsistent
 * DAG states.
 */
template<typename rlType>
void UserFunctionNode<rlType>::restoreMe( const RevBayesCore::DagNode* restorer )
{
    
    // We can no longer trust our value, so mark us as touched
    this->touched = true;
    
    // Dispatch call to downstream nodes
    this->restoreAffected();
}




template<class valueType>
void UserFunctionNode<valueType>::setMcmcMode(bool tf)
{
    
    // nothing to do
    
}


template<class valueType>
void UserFunctionNode<valueType>::setValueFromFile(const RevBayesCore::path& fn)
{
    
    //    *value = v;
    this->touch();
    
    throw RbException("Cannot set a user-function node from a file.");
}


template<class valueType>
void UserFunctionNode<valueType>::setValueFromString(const std::string &v)
{
    
    //    *value = v;
    this->touch();
    
    throw RbException("Cannot set a user-function node from a string.");
}


/**
 * Swap parent. We get this call just before an argument variable is being replaced by another variable.
 * We need to make sure we detach ourselves as a child of the old node and add ourselves as a child
 * to the new node. Actual replacement is taken care of by the variable wrapper of the variable.
 *
 * We need not worry about calls from cloneDAG, as we implement our own version of cloneDAG for this
 * DAG node class.
 */
template <typename rlType>
void UserFunctionNode<rlType>::swapParent(const RevBayesCore::DagNode* oldParent, const RevBayesCore::DagNode* newParent)
{
    // Find the old parent node
    const std::vector<const RevBayesCore::DagNode*>& parents = getParents();
    std::vector<const RevBayesCore::DagNode*>::const_iterator it;
    for ( it = parents.begin(); it != parents.end(); ++it )
    {
        if ( (*it) == oldParent )
        {
            // Detach old parent
            oldParent->removeChild( this );
            if ( oldParent->decrementReferenceCount() == 0 )
                delete oldParent;
            
            // Adopt new parent
            newParent->addChild( this );
            newParent->incrementReferenceCount();
            
            // Tell everybody we have been changed
            this->touch();
            
            // We are done
            return;
        }
    }
    
    // We did not find oldParent among our parents
    throw RbException( "Invalid attempt to swap non-parent" );
}


/**
 * Touch this node for recalculation. We only need to pass the message on
 * if we have not been touched before, conditional on all touched messages being
 * guaranteed to be followed either by a keep or a restore message so the DAG is
 * not in an inconsistent state. To be safe, we pass on the message regardless,
 * so that the touch propagates correctly regardless of the starting DAG state.
 */
template<typename rlType>
void UserFunctionNode<rlType>::touchMe( const RevBayesCore::DagNode* toucher, bool touchAll )
{
    
    // Touch myself
    this->touched = true;
    
    // Dispatch the touch message to downstream nodes
    this->touchAffected( touchAll );
}


/**
 * Update our value. We assume that the caller knows
 * that an update is really needed. We need to work
 * with pointers to be able to handle abstract value
 * types.
 */
template<typename rlType>
void UserFunctionNode<rlType>::update()
{
    // Update the return variable
    returnVariable = userFunction->executeCode();
    
    // We are clean!
    this->touched = false;
}


#endif
