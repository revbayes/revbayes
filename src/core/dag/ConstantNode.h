#ifndef ConstantNode_H
#define ConstantNode_H

#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "DagNodeMap.h"

namespace RevBayesCore {
    
    template<class valueType>
    class ConstantNode : public TypedDagNode<valueType> {
        
    public:
        ConstantNode(const std::string &n, valueType *v);
        ConstantNode(const ConstantNode<valueType> &c);                                                                                 //!< Copy constructor
        virtual                                            ~ConstantNode(void);                                                         //!< Virtual destructor
        
        void                                                bootstrap(void);                                                            //!< Bootstrap the current value of the node (applies only to stochastic nodes)
        virtual ConstantNode<valueType>*                    clone(void) const;                                                          //!< Create a clone of this node.
        DagNode*                                            cloneDAG(DagNodeMap &nodesMap, std::map<std::string, const DagNode* > &names) const; //!< Clone the entire DAG which is connected to this node
        void                                                getIntegratedParents(RbOrderedSet<DagNode *>& ip) const;
        double                                              getLnProbability(void);
        double                                              getLnProbabilityRatio(void);
        valueType&                                          getValue(void);
        const valueType&                                    getValue(void) const;
        bool                                                isConstant(void) const;                                                     //!< Is this DAG node constant?
        void                                                printStructureInfo(std::ostream &o, bool verbose=false) const;              //!< Print the structural information (e.g. name, value-type, distribution/function, children, parents, etc.)
        void                                                redraw(SimulationCondition c = SimulationCondition::MCMC);
        void                                                setMcmcMode(bool tf);                                                       //!< Set the modus of the DAG node to MCMC mode.
        void                                                setValue(const valueType &v);
        void                                                setValueFromFile(const path &dir);                                          //!< Set value from string.
        void                                                setValueFromString(const std::string &v);                                   //!< Set value from string.

    protected:
        void                                                getAffected(RbOrderedSet<DagNode *>& affected, const DagNode* affecter);    //!< Mark and get affected nodes
        void                                                keepMe(const DagNode* affecter);                                            //!< Keep value of this and affected nodes
        void                                                restoreMe(const DagNode *restorer);                                         //!< Restore value of this nodes
        void                                                touchMe(const DagNode *toucher, bool touchAll);                             //!< Tell affected nodes value is reset
        
    private:
        // members
        valueType*                                          value;
        
    };
    
    
}


#include "RbException.h"


template<class valueType>
RevBayesCore::ConstantNode<valueType>::ConstantNode(const std::string &n, valueType *v) : TypedDagNode<valueType>( n ),
    value( v )
{
    
    this->type = DagNode::CONSTANT;
    
}


template<class valueType>
RevBayesCore::ConstantNode<valueType>::ConstantNode(const ConstantNode<valueType> &c) : TypedDagNode<valueType>( c ),
    value( Cloner<valueType, IsDerivedFrom<valueType, Cloneable>::Is >::createClone( *c.value ) )
{
    
    this->type = DagNode::CONSTANT;
    
}

template<class valueType>
RevBayesCore::ConstantNode<valueType>::~ConstantNode( void )
{
    
    // we own the value so we need to delete it here
    delete value;
    
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::bootstrap( void )
{
    // nothing to do
}


/* Clone this node by creating an independent copy of the value. */
template<class valueType>
RevBayesCore::ConstantNode<valueType>* RevBayesCore::ConstantNode<valueType>::clone( void ) const
{
    
    return new ConstantNode<valueType>( *this );
    
}


/** Cloning the entire graph only involves children for a constant node */
template<class valueType>
RevBayesCore::DagNode* RevBayesCore::ConstantNode<valueType>::cloneDAG( DagNodeMap& newNodes, std::map<std::string, const DagNode* > &names ) const
{
    
    if ( newNodes.find( this ) != newNodes.end() )
    {
        return ( newNodes[ this ] );
    }
    
    // just for self checking purposes we keep track of the names for the variables we already cloned
    if ( this->name != "" )
    {
        // check if we already added a variable with this name
        std::map<std::string, const DagNode* >::const_iterator n = names.find( this->name );
        if ( n == names.end() )
        {
            // no, we haven't cloned a variable with this name before
            names[ this->name ] = this;
        }
        else
        {
            /*
            std::cerr << "Cloning a DAG node with name '" << this->name << "' again, doh!" << std::endl;
             */
        }
    }
    
    // Make pristine copy
    ConstantNode* copy = clone();
    newNodes[ this ] = copy;
    
    // Make sure the children clone themselves
    std::vector<DagNode*> children_to_clone = this->getChildren();
    for ( std::vector<DagNode* >::const_iterator i = children_to_clone.begin(); i != children_to_clone.end(); i++ )
    {
        DagNode *the_node = *i;
        std::string n = the_node->getName();
        the_node->cloneDAG( newNodes, names );
    }
    
    return copy;
}


/**
 * Get the affected nodes.
 * This call is started by the parent and since we don't have one this is a dummy implementation!
 */
template<class valueType>
void RevBayesCore::ConstantNode<valueType>::getAffected(RbOrderedSet<DagNode *> & /*affected*/, const DagNode* /*affecter*/)
{
    
    // do nothing
    throw RbException("You should never call getAffected() of a constant node!!!");
    
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::getIntegratedParents(RbOrderedSet<DagNode *>& ip) const
{
    // do nothing
}



template<class valueType>
double RevBayesCore::ConstantNode<valueType>::getLnProbability( void )
{
    
    return 0.0;
}


template<class valueType>
double RevBayesCore::ConstantNode<valueType>::getLnProbabilityRatio( void )
{
    
    return 0.0;
}


template<class valueType>
valueType& RevBayesCore::ConstantNode<valueType>::getValue( void )
{
    
    return *value;
}


template<class valueType>
const valueType& RevBayesCore::ConstantNode<valueType>::getValue( void ) const
{
    
    return *value;
}


template<class valueType>
bool RevBayesCore::ConstantNode<valueType>::isConstant( void ) const
{
    
    return true;
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::keepMe( const DagNode* /*affecter*/ )
{
    // nothing to do
}


template<class valueType>
/** Print struct for user */
void RevBayesCore::ConstantNode<valueType>::printStructureInfo(std::ostream &o, bool verbose) const
{
    
    o << "_dagType      = Constant node" << std::endl;
    o << "_children     = ";
    this->printChildren(o, 16, 70, verbose);
    o << std::endl;
    
    if ( verbose == true )
    {
        o << "_dagNode      = " << this->name << " <" << this << ">" << std::endl;
        o << "_refCount     = " << this->getReferenceCount() << std::endl;
    }
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::redraw( SimulationCondition c )
{
    // nothing to do
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::restoreMe( const DagNode * /*restorer*/ )
{
    // nothing to do
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::setMcmcMode(bool tf)
{
    
    // nothing to do
    
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::setValue(valueType const &v)
{
    
    *value = v;
    this->touch();
    
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::setValueFromFile(const RevBayesCore::path &dir)
{
    
    Serializer<valueType, IsDerivedFrom<valueType, RevBayesCore::Serializable>::Is >::ressurectFromFile( value, dir, this->getName() );
    this->touch();
    
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::setValueFromString(const std::string &v)
{
    Serializer<valueType, IsDerivedFrom<valueType, RevBayesCore::Serializable>::Is >::ressurectFromString( value, v );
    this->touch();
    
}


template<class valueType>
void RevBayesCore::ConstantNode<valueType>::touchMe( const DagNode * /*toucher*/, bool /*touchAll*/ )
{
    // nothing to do
}

#endif

