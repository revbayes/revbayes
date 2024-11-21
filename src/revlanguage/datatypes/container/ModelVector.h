#ifndef ModelVector_H
#define ModelVector_H

#include "RlContainer.h"
#include "ModelObject.h"
#include "RbVector.h"

#include <iostream>
#include <vector>

namespace RevLanguage {
    
    /**
     * @brief ModelVector: templated class for Rev vectors of model objects
     */
    template <typename rlType>
    class ModelVector : public ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >, public Container {
        
    public:
        typedef typename rlType::valueType                      elementType;
        typedef typename RevBayesCore::RbVector<elementType>    valueType;
        typedef typename valueType::iterator                    iterator;
        typedef typename valueType::const_iterator              const_iterator;
        
                                                                ModelVector(void);                                          //!< Default constructor
                                                                ModelVector(const valueType& v);                            //!< Constructor from vector of values
                                                                ModelVector(RevBayesCore::TypedDagNode<valueType>* n);      //!< Constructor from value node

        // STL-like vector functions provided here
        const elementType&                                      operator[](size_t index) const;                             //!< Subscript operator to internal value of Rev element, not allowing assignment
        void                                                    push_back(const rlType& x);                                 //!< Push Rev object element onto back
        void                                                    push_back(const RevObject &x);                              //!< Append element to end
        void                                                    push_back(const elementType& x);                            //!< Push internal value of Rev object element onto back
        
        // Basic utility functions you have to override
        virtual ModelVector<rlType>*                            clone(void) const;                                          //!< Clone object
        static const std::string&                               getClassType(void);                                         //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                     //!< Get class type spec
        virtual const TypeSpec&                                 getTypeSpec(void) const;                                    //!< Get the object type spec of the instance
        
 
        // Type conversion functions
        RevObject*                                              convertTo(const TypeSpec& type) const;                      //!< Convert to requested type
        virtual double                                          isConvertibleTo(const TypeSpec& type, bool once) const;     //!< Is this object convertible to the requested type?

        // Member object functions
        virtual RevPtr<RevVariable>                             executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found); //!< Map member methods to internal methods
        
        // Container functions provided here
        virtual rlType*                                         getElement(size_t idx) const;                               //!< Get element variable (single index)

        // ModelVector functions: override if you do not want to support these in-place algorithms
        virtual void                                            clear(void);                                                //!< Clear the vector
        virtual size_t                                          size(void) const;                                           //!< Size of the vector
        virtual void                                            sort(void);                                                 //!< Sort vector
        virtual void                                            unique(void);                                               //!< Remove consecutive duplicates
        void                                                    printValue(std::ostream& o, bool user) const;               //!< Print value for user

        std::string                                             getGuiName(void) { return "Vector"; }
        std::string                                             getGuiUnicodeSymbol(void) { return "v"; }
        std::string                                             getGuiInfo(void) { return ""; }

    private:
        void                                                    initMethods(void);
        
        struct comparator {
            bool operator() (elementType & A, elementType & B) const { return ( A < B ); }
        } myComparator;
    
    };
}


#include "ArgumentRule.h"
#include "DeterministicNode.h"
#include "MethodTable.h"
#include "Natural.h"
#include "Probability.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "TypeSpec.h"
#include "Workspace.h"
#include "WorkspaceVector.h"

#include <algorithm>

using namespace RevLanguage;


/**
 * Default constructor: construct an empty vector. We do not need to store the element type
 * because this can be accessed from our Rev type specification.
 */
template <typename rlType>
ModelVector<rlType>::ModelVector( void ) :
ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >( new RevBayesCore::RbVector<typename rlType::valueType>() )
{
    
    initMethods();
    
}


/**
 * Constructor from a valueType object, that is, from an STL vector containing instances
 * of the internal value type of the Rev object elements of the model vector.
 */
template <typename rlType>
ModelVector<rlType>::ModelVector( const valueType &v ) :
    ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >( v.clone() )
{
    
    initMethods();
    
}


/**
 * Constructor from a DAG node appropriate for the model vector
 * object.
 */
template <typename rlType>
ModelVector<rlType>::ModelVector( RevBayesCore::TypedDagNode<valueType> *n ) :
    ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >( n )
{
    
    initMethods();
    
}


/**
 * Subscript operator, provided for convenience. Note that we should
 * not allow the user to assign to our elements, so we need to give
 * out a copy and not a reference to the element. This makes it
 * possible to implement this as a const function returning a non-const
 * variable. The return of a copy is ensured by not specifying a reference
 * (&) on the return type.
 */
template <typename rlType>
const typename rlType::valueType& ModelVector<rlType>::operator[]( size_t index ) const
{
    return this->getValue()[ index ];
}


/**
 * Clear the vector.
 */
template <typename rlType>
void ModelVector<rlType>::clear( void )
{
    this->dag_node->getValue().clear();
}


/** Get a type-safe clone of the object */
template <typename rlType>
ModelVector<rlType>* ModelVector<rlType>::clone() const
{
    return new ModelVector<rlType>( *this );
}


/**
 * Convert to object of another type. Here we use the setElements function
 * of the Container base class to do generic type conversion in all cases
 * where the elements are individually convertible to the desired element
 * type. This is not done automatically for us because of the templating.
 * A vector of RealPos, for example, does not inherit from a vector of Real,
 * which means that a vector of RealPos is not a specialized vector of Real
 * in the C++ sense (or in the Rev sense).
 */
template <typename rlType>
RevObject* ModelVector<rlType>::convertTo(const TypeSpec &type) const
{
    
    // First check that we are not asked to convert to our own type
    if ( type == getClassTypeSpec() )
    {
        return this->clone();
    }
    
    // Test whether we want to convert to another generic model vector
    if ( type.getParentType() == getClassTypeSpec().getParentType() )
    {
        // We are both model vectors. Rely on generic code to cover all allowed conversions

        // First generate an empty model vector of the desired type
        RevObject *emptyContainer = Workspace::userWorkspace().makeNewDefaultObject( type.getType() );
        Container *theConvertedContainer = dynamic_cast<Container*>( emptyContainer );
        
        // test if the cast succeeded
        if ( theConvertedContainer == NULL )
        {
            throw RbException() << "Could not convert a container of type " << this->getClassType() << " to a container of type " <<  type.getType() ;
        }

        for ( typename RevBayesCore::RbConstIterator<elementType> i = this->getValue().begin(); i != this->getValue().end(); ++i )
        {
            
            rlType orgElement = rlType( *i );
            if ( type.getElementTypeSpec() != NULL  && orgElement.isType( *type.getElementTypeSpec() ) == true )
            {
                theConvertedContainer->push_back( orgElement );
            }
            else
            {
                RevObject *convObj = orgElement.convertTo( *type.getElementTypeSpec() );
                theConvertedContainer->push_back( *convObj );
                delete convObj;
            }
            
        }
        
        // Now return the converted container object
        return emptyContainer;
    }
    else if ( type == WorkspaceVector<RevObject>::getClassTypeSpec() )
    {
        // create an empty container
        WorkspaceVector<RevObject> *theConvertedContainer = new WorkspaceVector<RevObject>();

        if ( this->getDagNode()->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC )
        {
            
            std::vector<const RevBayesCore::DagNode*> args = this->getDagNode()->getParents();
            
            for ( std::vector<const RevBayesCore::DagNode*>::iterator i = args.begin(); i != args.end(); ++i )
            {
                RevBayesCore::DagNode* node = const_cast<RevBayesCore::DagNode*>(*i);
                RevBayesCore::TypedDagNode<elementType>* tnode = static_cast<RevBayesCore::TypedDagNode<elementType>* >( node );
                rlType orgElement = rlType( tnode );
                theConvertedContainer->push_back( orgElement );
                
            }
            
        }
        else
        {
            
            for ( typename RevBayesCore::RbConstIterator<elementType> i = this->getValue().begin(); i != this->getValue().end(); ++i )
            {
            
                rlType orgElement = rlType( *i );
                theConvertedContainer->push_back( orgElement );
            
            }
            
        }
        
        return theConvertedContainer;
    }
    
    
    // Call the base class if all else fails. This will eventually throw an error if the type conversion is not supported.
    return this->ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >::convertTo( type );
}


/**
 * Map calls to member methods.
 */
template <typename rlType>
RevPtr<RevVariable> ModelVector<rlType>::executeMethod( std::string const &name, const std::vector<Argument> &args, bool &found )
{
    
    if ( name == "append" )
    {
        found = true;
        
        // Check whether the DAG node is actually a constant node
        if ( this->dag_node->isConstant() == false )
        {
            throw RbException( "Only constant variables can be appended." );
        }
        
        RevBayesCore::RbVector<typename rlType::valueType> &v = this->dag_node->getValue();
        
        if ( args[0].getVariable()->getRevObject().isType( ModelVector<rlType>::getClassTypeSpec() ) )
        {
            const ModelVector<rlType> &v_x = static_cast<const ModelVector<rlType>&>( args[0].getVariable()->getRevObject() );
            const RevBayesCore::RbVector<typename rlType::valueType> &x = v_x.getValue();
            for (size_t i = 0; i < x.size(); ++i )
            {
                v.push_back( x[i] );
            }
        }
        else
        {
            const rlType &rl_x = static_cast<const rlType&>( args[0].getVariable()->getRevObject() );
            const typename rlType::valueType &x = rl_x.getValue();
            v.push_back( x );
        }
        
        return NULL;
    }
    else if ( name == "contains" )
    {
        found = true;
        
        const rlType &rl_x = static_cast<const rlType&>( args[0].getVariable()->getRevObject() );
        const typename rlType::valueType &x = rl_x.getValue();
        const RevBayesCore::RbVector<typename rlType::valueType> &v = this->dag_node->getValue();
        for (size_t i = 0; i < v.size(); ++i )
        {
            if ( v[i] == x )
            {
                return RevPtr<RevVariable>( new RevVariable( new RlBoolean( true ), "" ) );
            }
        }
        
        // return a new RevVariable with the size of this container
        return RevPtr<RevVariable>( new RevVariable( new RlBoolean( false ), "" ) );
    }
    else if ( name == "erase" )
    {
        found = true;
        
        // Check whether the DAG node is actually a constant node
        if ( this->dag_node->isConstant() == false )
        {
            throw RbException( "We can only erase from constant variables." );
        }
        
        RevBayesCore::RbVector<typename rlType::valueType> &v = this->dag_node->getValue();
        
        if ( args[0].getVariable()->getRevObject().isType( ModelVector<rlType>::getClassTypeSpec() ) )
        {
            const ModelVector<rlType> &v_x = static_cast<const ModelVector<rlType>&>( args[0].getVariable()->getRevObject() );
            const RevBayesCore::RbVector<typename rlType::valueType> &x = v_x.getValue();
            for (size_t i = 0; i < x.size(); ++i )
            {
                RevBayesCore::RbIterator<typename rlType::valueType> pos = v.find( x[i] );
                if ( pos != v.end() )
                {
                    v.erase( pos );
                }
            }
        }
        else
        {
            const rlType &rl_x = static_cast<const rlType&>( args[0].getVariable()->getRevObject() );
            const typename rlType::valueType &x = rl_x.getValue();
            RevBayesCore::RbIterator<typename rlType::valueType> pos = v.find( x );
            if ( pos != v.end() )
            {
                v.erase( pos );
            }
        }
        
        return NULL;
    }
    else if ( name == "find" )
    {
        found = true;
        
        const rlType &rl_x = static_cast<const rlType&>( args[0].getVariable()->getRevObject() );
        const typename rlType::valueType &x = rl_x.getValue();
        const RevBayesCore::RbVector<typename rlType::valueType> &v = this->dag_node->getValue();
        for (size_t i = 0; i < v.size(); ++i )
        {
            if ( v[i] == x )
            {
                return RevPtr<RevVariable>( new RevVariable( new Natural( i+1 ), "" ) );
            }
        }
        
        // return a new RevVariable with the size of this container
        return RevPtr<RevVariable>( new RevVariable( new Natural( -1 ), "" ) );
    }
    else if ( name == "size" )
    {
        found = true;
        
        // return a new RevVariable with the size of this container
        return RevPtr<RevVariable>( new RevVariable( new Natural( size() ), "" ) );
    }
    else if ( name == "sort" )
    {
        found = true;
        
        // Check whether the DAG node is actually a constant node
        if ( this->dag_node->isConstant() == false )
        {
            throw RbException( "Only constant variables can be sorted." );
        }
        sort();
        
        return NULL;
    }
    else if ( name == "unique" )
    {
        found = true;
        
        // Check whether the DAG node is actually a constant node
        if ( this->dag_node->isConstant() == false )
        {
            throw RbException( "Only constant variables can be made unique." );
        }
        unique();
        
        return NULL;
    }
    else if ( name == "[]" )
    {
        found = true;
        
        long index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        return RevPtr<RevVariable>( new RevVariable( getElement( index ) ) );
    }
    
    return ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >::executeMethod( name, args, found );
}


/**
 * Get Rev type of object. This is the Rev object element type followed by
 * a single set of square brackets. This provides a nice and convenient way
 * of specifying generic types of vectors for all Rev object types.
 */
template <typename rlType>
const std::string& ModelVector<rlType>::getClassType(void)
{
    static std::string rev_type = rlType::getClassType() + "[]";
    
	return rev_type;
}


/**
 * Get class type spec describing type of object. Note that we have to
 * use the special version of the TypeSpec constructor because we derive
 * directly from ModelContainer and have the special "rlType[]" type
 * specification.
 */
template <typename rlType>
const RevLanguage::TypeSpec& ModelVector<rlType>::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >::getClassTypeSpec(), &rlType::getClassTypeSpec() );

	return rev_type_spec;
}


template <typename rlType>
rlType* ModelVector<rlType>::getElement(size_t idx) const
{
    return new rlType( this->getValue()[ idx ] );
}


/** Get the type spec (dynamic version) */
template <typename rlType>
const TypeSpec& ModelVector<rlType>::getTypeSpec(void) const
{
    return getClassTypeSpec();
}


/**
 * Initialize the methods.
 */
template <typename rlType>
void ModelVector<rlType>::initMethods( void )
{

    ArgumentRules* append_arg_rules = new ArgumentRules();
    
    std::vector<TypeSpec> appendValueTypes;
    appendValueTypes.push_back( rlType::getClassTypeSpec() );
    appendValueTypes.push_back( ModelVector<rlType>::getClassTypeSpec() );
    append_arg_rules->push_back( new ArgumentRule( "x", appendValueTypes, "The element.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "append", RlUtils::Void, append_arg_rules ) );

    ArgumentRules* contains_arg_rules = new ArgumentRules();
    contains_arg_rules->push_back( new ArgumentRule( "x", rlType::getClassTypeSpec(), "The element.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "contains", RlBoolean::getClassTypeSpec(), contains_arg_rules ) );

    ArgumentRules* erase_arg_rules = new ArgumentRules();
    std::vector<TypeSpec> erase_value_types;
    erase_value_types.push_back( rlType::getClassTypeSpec() );
    erase_value_types.push_back( ModelVector<rlType>::getClassTypeSpec() );
    erase_arg_rules->push_back( new ArgumentRule( "x", erase_value_types, "The element that you want to erase.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "erase", RlUtils::Void, erase_arg_rules) );

    ArgumentRules* find_arg_rules = new ArgumentRules();
    find_arg_rules->push_back( new ArgumentRule( "x", rlType::getClassTypeSpec(), "The element.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "find", RlBoolean::getClassTypeSpec(), find_arg_rules ) );
    
    ArgumentRules* sizeArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );
    
    ArgumentRules* sortArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "sort", RlUtils::Void, sortArgRules) );
    
    ArgumentRules* uniqueArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "unique", RlUtils::Void, uniqueArgRules) );
    
    ArgumentRules* elementArgRules = new ArgumentRules();
    elementArgRules->push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The index of the element.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", rlType::getClassTypeSpec(), elementArgRules ) );
    
}


/**
 * In this function we check whether this type is convertible to some other
 * Rev object type. Here we focus entirely on supporting conversion to
 * other generic vectors with compatible elements. This is not done automatically
 * because of the templating: a vector of RealPos does not inherit from a vector
 * of Real, for example.
 */
template <typename rlType>
double ModelVector<rlType>::isConvertibleTo( const TypeSpec& type, bool once ) const
{

    if ( once == true && type.getParentType() == getClassTypeSpec().getParentType() )
    {
        // We want to convert to another model vector
        if ( getClassType() == "RealPos[]" && type.getType() == "Real[]" )
        {
            return this->getValue().size() * 0.2;
        }
        if ( getClassType() == "Natural[]" && type.getType() == "Integer[]" )
        {
            return this->getValue().size() * 0.2;
        }

        // Simply check whether our elements can be converted to the desired element type
        typename RevBayesCore::RbConstIterator<elementType> i;
        double penalty = 0.0;
        double sum = 0.0;
        for ( i = this->getValue().begin(); i != this->getValue().end(); ++i )
        {
            const elementType& org_internal_element = *i;
            rlType org_element = rlType( org_internal_element );

            // Test whether this element is already of the desired element type or can be converted to it
            if ( type.getElementTypeSpec() == NULL )
            {
                // we cannot convert this element
                return -1;
            }
            
            if ( org_element.getTypeSpec() != *type.getElementTypeSpec() )
            {
            
                double element_penalty = org_element.isConvertibleTo( *type.getElementTypeSpec(), once );
                if ( element_penalty == -1 )
                {
                    // we cannot convert this element
                    return -1;
                }
                penalty += element_penalty;
            }

            if( type.getType() == "Simplex" )
            {
                Probability* ro = static_cast<Probability*>( org_element.convertTo( Probability::getClassTypeSpec() ) );
                sum += ro->getValue();
                delete ro;
            }
        }

        // we cannot convert to a Simplex if the elements are not normalized
        if ( type.getType() == "Simplex" && RevBayesCore::RbMath::compEssentiallyEqual(sum, 1.0) == false )
        {
            return -1;
        }

        return penalty;
    }
    else if ( type == WorkspaceVector<RevObject>::getClassTypeSpec() )
    {
        // yes we can
        return 0.0;
    }

    return ModelObject<RevBayesCore::RbVector<typename rlType::valueType> >::isConvertibleTo( type, once );
}


/**
 * Add an element to the end of the vector.
 */
template <typename rlType>
void ModelVector<rlType>::push_back(const elementType &x)
{
    return this->dag_node->getValue().push_back( x );
}


/**
 * Add an element to the end of the vector.
 */
template <typename rlType>
void ModelVector<rlType>::push_back(const rlType &x)
{
    return this->dag_node->getValue().push_back( x.getValue() );
}


/**
 * Push a Rev object element onto the back of the vector.
 */
template <typename rlType>
void ModelVector<rlType>::push_back( const RevObject &x )
{
    
    // cast the object
    const rlType *x_converted = dynamic_cast< const rlType* >( &x );
    
    if ( x_converted == NULL )
    {
        throw RbException() << "Could not append an element of type " << x.getType() << " to a vector of type " <<  this->getType() ;
    }
    
    // Push it onto the back of the elements vector
    this->push_back( *x_converted );
}


/**
 * Print the value of the vector, respecting the formatting of the
 * model object elements. We do this by retrieving the elements, one
 * by one, and printing them using their own printValue implementation.
 * Among other things, this takes care of proper formatting.
 *
 * We make a perfectly safe const cast here, since we only utilize the
 * const printValue function of the element.
 */
template <typename rlType>
void ModelVector<rlType>::printValue( std::ostream& o, bool user ) const
{
    
    this->getDagNode()->printValue( o, ",", -1, true, user, true );

}


/**
 * Size of the vector.
 */
template <typename rlType>
size_t ModelVector<rlType>::size( void ) const
{
    return this->dag_node->getValue().size();
}


/**
 * Sort the vector.
 */
template <typename rlType>
void ModelVector<rlType>::sort( void )
{
    this->dag_node->getValue().sort();
}

                                       
/**
 * Remove duplicates and resize the vector. This is a mutable member function,
 * meaning that it modifies this instance.
 *
 */
template <typename rlType>
void ModelVector<rlType>::unique(void)
{
    
    sort();
    valueType uniqueVector;
    uniqueVector.push_back (this->dag_node->getValue()[0]);
    for (size_t i = 1 ; i<this->size() ; i++)
    {
        if ( this->dag_node->getValue()[i] != this->dag_node->getValue()[i-1] )
        {
            uniqueVector.push_back(this->dag_node->getValue()[i]);
        }
    }
    
    this->clear();

    this->dag_node->getValue() = uniqueVector;
    
}


#endif
