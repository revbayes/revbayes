#ifndef ModelObject_H
#define ModelObject_H

#include "AbstractModelObject.h"
#include "MethodTable.h"
#include "RevObject.h"
#include "TypedDagNode.h"
#include "UserFunction.h"

namespace RevLanguage {
    
    template <typename rbType>
    class ModelObject : public AbstractModelObject {
    
    public:
        virtual                                ~ModelObject(void);                                                          //!< Virtual destructor

        ModelObject&                            operator=(const ModelObject& x);                                            //!< Assignment operator
        
        // The value type definition
        typedef rbType                          valueType;
       
        // Basic utility functions you have to override
        virtual ModelObject*                    clone(void) const = 0;                                                      //!< Clone object
        static const std::string&               getClassType(void);                                                         //!< Get Rev type (static)
        static const TypeSpec&                  getClassTypeSpec(void);                                                     //!< Get Rev type spec (static)
        virtual const TypeSpec&                 getTypeSpec(void) const = 0;                                                //!< Get Rev type spec (instance)
    
        // Utility functions you might want to override
        virtual RevPtr<RevVariable>             executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found);  //!< Execute member functions
        
        // Basic utility functions you should not have to override
        bool                                    isAssignable(void) const;                                                   //!< Is object or upstream members assignable?
        bool                                    isConstant(void) const;                                                     //!< Is this variable and the internally stored deterministic node constant?
        void                                    makeConstantValue(void);                                                    //!< Convert to constant object
        void                                    makeConversionValue(RevPtr<RevVariable> var);                                  //!< Convert to conversion object
        ModelObject<rbType>*                    makeIndirectReference(void);                                                //!< Make reference to object
        void                                    makeUserFunctionValue(UserFunction* fxn);                                   //!< Convert to user-defined Rev function object
        void                                    setDagNode(RevBayesCore::DagNode *newNode);                                 //!< Set or replace the internal dag node (and keep me)
        void                                    setName(const std::string &n);                                              //!< Set the name of the variable (if applicable)
        void                                    replaceVariable(RevObject *newVar);                                         //!< Replace the internal DAG node (and prepare to replace me...)
        
        // Getters and setters
        RevBayesCore::TypedDagNode<rbType>*     getDagNode(void) const;                                                     //!< Get the internal DAG node
        bool                                    hasDagNode(void) const;                                                     //!< Get the internal DAG node
        virtual const rbType&                   getValue(void) const;                                                       //!< Get the value (const)
        virtual rbType&                         getValue(void);                                                             //!< Get the value (non-const)
        void                                    setValue(rbType *x);                                                        //!< Set new constant value
        virtual void                            printValue(std::ostream& o, bool user) const;                               //!< Print value
        virtual void                            printValue(std::ostream& o) const { printValue(o, true); };                 //!< Print value overload

        virtual std::string                     getGuiName(void) { return ""; }
        virtual std::string                     getGuiUnicodeSymbol(void) { return ""; }
        virtual std::string                     getGuiInfo(void) { return ""; }
        
    protected:
        ModelObject(void);
        ModelObject(rbType *v);
        ModelObject(RevBayesCore::TypedDagNode<rbType> *v);
        ModelObject(const ModelObject &v);
        
        RevBayesCore::TypedDagNode<rbType>*     dag_node;

    };
    
}


#include "AbstractCharacterData.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Cloner.h"
#include "ConstantNode.h"
#include "ContinuousCharacterData.h"
#include "Func__conversion.h"
#include "IndirectReferenceFunction.h"
#include "MemberProcedure.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlUtils.h"
#include "StochasticNode.h"
#include "UserFunctionNode.h"
#include "RevVariable.h"
#include "Workspace.h"

#include <cassert>
#include <cmath>

template <typename rbType>
RevLanguage::ModelObject<rbType>::ModelObject() : AbstractModelObject(),
    dag_node( NULL )
{
}



template <typename rbType>
RevLanguage::ModelObject<rbType>::ModelObject(rbType *v) : AbstractModelObject(),
    dag_node( new ConstantNode<rbType>("",v) )
{
    // increment the reference count to the value
    dag_node->incrementReferenceCount();
    
}



template <typename rbType>
RevLanguage::ModelObject<rbType>::ModelObject(RevBayesCore::TypedDagNode<rbType> *v) : AbstractModelObject(),
    dag_node( v )
{
    // increment the reference count to the value
    dag_node->incrementReferenceCount();
    
    // add the DAG node member methods
    // note that this is safe case because all DAG nodes are member objects
    const MethodTable &dagMethods = dynamic_cast<RevMemberObject*>( dag_node )->getMethods();
    methods.insertInheritedMethods( dagMethods );

}



template <typename rbType>
RevLanguage::ModelObject<rbType>::ModelObject(const ModelObject &v) :
    AbstractModelObject( v ),
    dag_node( NULL )
{
    
    if ( v.dag_node != NULL )
    {
        
        dag_node = v.dag_node->clone();
        
        // increment the reference count to the value
        dag_node->incrementReferenceCount();
    
    }
    
}


template <typename rbType>
RevLanguage::ModelObject<rbType>::~ModelObject() 
{
    
    // free the old value
    if ( dag_node != NULL )
    {
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
    }
    
}


template <typename rbType>
RevLanguage::ModelObject<rbType>& RevLanguage::ModelObject<rbType>::operator=(const ModelObject &v)
{
    
    if ( this != &v ) 
    {
        // delegate to base class
        AbstractModelObject::operator=( v );
        
        // free the old value
        if ( dag_node != NULL )
        {
            if ( dag_node->decrementReferenceCount() == 0 )
            {
                delete dag_node;
            }
            
            dag_node = NULL;
        }
        
        // create own copy
        if ( v.dag_node != NULL )
        {
            dag_node = v.dag_node->clone();
            
            // increment the reference count to the value
            dag_node->incrementReferenceCount();
        
        }
    }
    
    return *this;
}


/* Map calls to member methods */
template <typename rbType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::ModelObject<rbType>::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if ( name == "methods" )
    {
        found = true;
        
        // just print the method names (including inherited methods)
        const MethodTable &m = getMethods();
        m.printValue(std::cout, true);
        
        return NULL;
    }
    
    RevMemberObject * rmo = dynamic_cast<RevMemberObject *>( dag_node );
    
    if ( rmo != NULL )
    {
        RevPtr<RevVariable> retVal = rmo->executeMethod(name, args, found);

        if ( found == true)
        {
            return retVal;
        }
    }
    
    return RevObject::executeMethod( name, args, found );
}


template <typename rbType>
const std::string& RevLanguage::ModelObject<rbType>::getClassType(void)
{
    
    static std::string rev_type = "ModelObject";
    
	return rev_type;
}


/** Get class type spec describing type of object */
template <typename rbType>
const RevLanguage::TypeSpec& RevLanguage::ModelObject<rbType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &RevObject::getClassTypeSpec() );
    
	return rev_type_spec;
}



template <typename rbType>
RevBayesCore::TypedDagNode<rbType>* RevLanguage::ModelObject<rbType>::getDagNode( void ) const
{
    return dag_node;
}

template <typename rbType>
bool RevLanguage::ModelObject<rbType>::hasDagNode( void ) const
{
    return (dag_node != nullptr);
}

template <typename rbType>
const rbType& RevLanguage::ModelObject<rbType>::getValue( void ) const
{
    
    if ( dag_node == NULL )
    {
        throw RbException( "Invalid attempt to get value from an object with NULL DAG node" );
    }
    
    return dag_node->getValue();
}


template <typename rbType>
rbType& RevLanguage::ModelObject<rbType>::getValue( void )
{
    
    if ( dag_node == NULL )
        throw RbException( "Invalid attempt to get value from an object with NULL DAG node" );
    
    return dag_node->getValue();
}


/**
 * Is the object or any of its upstream members or elements
 * modifiable by the user through assignment? We simply ask
 * our DAG node.
 */
template <typename rbType>
bool RevLanguage::ModelObject<rbType>::isAssignable( void ) const
{
    if ( dag_node == NULL )
        return false;
    
    return dag_node->isAssignable();
}


template <typename rbType>
bool RevLanguage::ModelObject<rbType>::isConstant( void ) const
{
    
    return dag_node->isConstant();
}


template <typename rbType>
void RevLanguage::ModelObject<rbType>::makeConstantValue( void )
{
    
    if ( dag_node == NULL )
    {
        throw RbException("Cannot convert a variable without value to a constant value.");
    }
    else if ( dag_node->getDagNodeType() != RevBayesCore::DagNode::CONSTANT )
    {
        RevBayesCore::ConstantNode<rbType>* new_node = new ConstantNode<rbType>(dag_node->getName(), RevBayesCore::Cloner<rbType, IsDerivedFrom<rbType, RevBayesCore::Cloneable>::Is >::createClone( dag_node->getValue() ) );
        dag_node->replace(new_node);
        
        // delete the value if there are no other references to it.
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
        
        dag_node = new_node;
        
        // increment the reference counter
        dag_node->incrementReferenceCount();
    }
    
}


/**
 * Make an indirect reference to the variable. This is appropriate for the contexts
 * where the object occurs on the righ-hand side of expressions like a := b
 */
template <typename rbType>
RevLanguage::ModelObject<rbType>* RevLanguage::ModelObject<rbType>::makeIndirectReference(void)
{
//    IndirectReferenceNode< RevLanguage::ModelObject<rbType> >* newNode = new IndirectReferenceNode< RevLanguage::ModelObject<rbType> >( "", this->getDagNode() );
    
    RevBayesCore::IndirectReferenceFunction< rbType > *func = new RevBayesCore::IndirectReferenceFunction<rbType>( this->getDagNode() );
    RevBayesCore::DeterministicNode< rbType >* newNode = new RevBayesCore::DeterministicNode< rbType >( "", func );
    
    RevLanguage::ModelObject< rbType >* newObj = this->clone();
    
    const std::vector<RevBayesCore::Move*>& mvs = newObj->getDagNode()->getMoves();
    while ( mvs.empty() == false )
    {
        newObj->getDagNode()->removeMove( *mvs.begin() );
    }
    
    newObj->setDagNode( newNode );
    
    return newObj;
}


/** Convert a model object to a deterministic object, the value of which is determined by a user-defined Rev function */
template <typename rbType>
void RevLanguage::ModelObject<rbType>::makeUserFunctionValue( UserFunction* fxn )
{
    UserFunctionNode< ModelObject<rbType> >*  detNode = new UserFunctionNode< ModelObject<rbType> >( "", fxn );
    
    // Signal replacement and delete the value if there are no other references to it.
    if ( dag_node != NULL )
    {
        dag_node->replace( detNode );
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
        
    }
    
    // Shift the actual node
    dag_node = detNode;
    
    // Increment the reference counter
    dag_node->incrementReferenceCount();
}


/**
 * Print value for user. The DAG node pointer may be NULL, in which
 * case we print "NA".
 */
template <typename rbType>
void RevLanguage::ModelObject<rbType>::printValue(std::ostream &o, bool user) const
{
    if ( dag_node == NULL )
    {
        o << "NA";
    }
    else
    {
        dag_node->printValue( o, "", -1, true, user, true );
    }
    
}


/** Copy name of variable onto DAG node, if it is not NULL */
template <typename rbType>
void RevLanguage::ModelObject<rbType>::setName(std::string const &n)
{
    if ( dag_node != NULL )
    {
        dag_node->setName( n );
    }
    
}


/**
 * Set dag node. We also accommodate the possibility of setting the DAG node to null.
 */
template <typename rbType>
void RevLanguage::ModelObject<rbType>::setDagNode(RevBayesCore::DagNode* newNode)
{
    
    // Take care of the old value node
    if ( dag_node != NULL )
    {
        if ( newNode != NULL )
        {
            newNode->setName( dag_node->getName() );
        }
        
        dag_node->replace(newNode);
        
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
        
    }
    
    // Set the new value node
    dag_node = static_cast< RevBayesCore::TypedDagNode<rbType>* >( newNode );
    
    // Increment the reference count to the new value node
    if ( dag_node != NULL )
    {
        dag_node->incrementReferenceCount();
    }
    
}


template <typename rbType>
void RevLanguage::ModelObject<rbType>::setValue(rbType *x)
{
    
    RevBayesCore::ConstantNode<rbType>* newNode;
    
    if ( dag_node == NULL )
    {
        newNode = new ConstantNode<rbType>("",x);
    }
    else
    {
        newNode = new ConstantNode<rbType>(dag_node->getName(),x);
        dag_node->replace(newNode);
        
        if ( dag_node->decrementReferenceCount() == 0 )
        {
            delete dag_node;
        }
        
    }
    
    dag_node = newNode;
    
    // increment the reference count to the value
    dag_node->incrementReferenceCount();
    
}


#endif
