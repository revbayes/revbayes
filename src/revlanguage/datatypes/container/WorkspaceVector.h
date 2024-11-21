#ifndef WorkspaceVector_H
#define WorkspaceVector_H

#include "RlContainer.h"
#include "RbConstIterator.h"
#include "RbVector.h"
#include "WorkspaceToCoreWrapperObject.h"

#include <iostream>
#include <vector>

namespace RevLanguage {
    
    /**
     * @brief WorkspaceVector: templated vector of workspace object elements
     *
     * WorkspaceVector is a container templated on the
     * Rev type of the elements of the vector, which are expected to be workspace objects
     * (inheriting from WorkspaceObject).
     *
     * The class gives access to an RbVector of clones of the internal value objects through
     * the getVectorRbPointer() method. The workspace vector can be set using either an
     * RbVector object or a vector of RevObject element pointers (vectorRlPtr). In the latter
     * case, we assume ownership of the elements; in the former, we generate independent clones
     * of the elements.
     *
     * For Rev objects that are not model objects or workspace objects, see RevObjectVector.
     */
    template <typename rlType>
    class WorkspaceVector : public WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >, public Container {
        
    public:

                                                    WorkspaceVector(void);                                              //!< Default constructor
                                                    WorkspaceVector(const RevBayesCore::RbVector<rlType>& v);           //!< Constructor from rb core vector
        
        // STL-like vector functions provided here
        rlType&                                     operator[](size_t index);                                           //!< Subscript operator
        const rlType&                               operator[](size_t index) const;                                     //!< Subscript operator (const)
        RevBayesCore::RbConstIterator<rlType>       end(void) const { /* NEEDS IMPLEMENATION!!! */ return RevBayesCore::RbConstIterator<rlType>(); }                                                    //!< Iterator to end
        void                                        erase(RevBayesCore::RbConstIterator<rlType> p) { /* NEEDS IMPLEMENATION!!! */ }                     //!< Erase an element
        RevBayesCore::RbConstIterator<rlType>       find(const rlType &x) const { /* NEEDS IMPLEMENATION!!! */ return RevBayesCore::RbConstIterator<rlType>(); }                                        //!< Find an element
        rlType&                                     pop_back(void);                                                     //!< Drop element at back
        rlType&                                     pop_front(void);                                                    //!< Drop element from front
        void                                        push_back(const RevObject &x);                                      //!< Append element to front
        void                                        push_front(const rlType &x);                                        //!< Append element to end
        size_t                                      size(void) const;                                                   //!< The size of the vector

        // Basic utility functions you have to override
        virtual WorkspaceVector<rlType>*            clone(void) const;                                                  //!< Clone object
        static const std::string&                   getClassType(void);                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                             //!< Get class type spec
        virtual const TypeSpec&                     getTypeSpec(void) const;                                            //!< Get the object type spec of the instance
        virtual RevPtr<RevVariable>                 executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found); //!< Map member methods to internal methods

        // Type conversion functions
        RevObject*                                  convertTo(const TypeSpec& type) const;                              //!< Convert to requested type
        virtual double                              isConvertibleTo(const TypeSpec& type, bool once) const;             //!< Is this object convertible to the requested type?

        // Container functions provided here
        virtual rlType*                             getElement(size_t idx) const;                                       //!< Get element variable (single index)

        // WorkspaceVector functions
        void                                        printValue(std::ostream& o, bool) const;                       //!< Print value for user

        
    private:
        void                                        initMethods(void);  //!< Initialize member methods.

    };

}


#include "Natural.h"
#include <algorithm>

using namespace RevLanguage;


/**
 * Default constructor. Construct empty vector.
 * We do not have to set the element type because that is
 * taken care of in getClassType.
 */
template <typename rlType>
WorkspaceVector<rlType>::WorkspaceVector( void ) :
    WorkspaceToCoreWrapperObject<RevBayesCore::RbVector< rlType > >()
{
    this->value = new RevBayesCore::RbVector< rlType >();
    
    initMethods();
}


/**
 * Construct from vector of pointers to internal value
 * objects. We assume here that we are responsible for
 * managing the associated memory. Here we need to make
 * copies of the original objects, so we delete those
 * after we are done.
 */
template <typename rlType>
WorkspaceVector<rlType>::WorkspaceVector( const RevBayesCore::RbVector<rlType>& v ) : WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >()
{
    // just to be safe
    delete this->value;
    
    this->value = v.clone();
    initMethods();
}


/**
 * Subscript operator, provided for convenience. Note that
 * there is no problem to give out non-const references
 * to our Rev object elements.
 */
template <typename rlType>
rlType& WorkspaceVector<rlType>::operator[]( size_t index )
{
    return (*this->value)[ index ];
}


/**
 * Subscript operator, provided for convenience. This is
 * the const version.
 */
template <typename rlType>
const rlType& WorkspaceVector<rlType>::operator[]( size_t index ) const
{
    return (*this->value)[ index ];
}


/** Return a type-safe clone of the object */
template <typename rlType>
WorkspaceVector<rlType>* WorkspaceVector<rlType>::clone() const
{
    return new WorkspaceVector<rlType>( *this );
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
RevObject* WorkspaceVector<rlType>::convertTo(const TypeSpec &type) const
{
    
    // First check that we are not asked to convert to our own type
    if ( type == getClassTypeSpec() )
    {
        return this->clone();
    }
    
    if ( type.getParentTypeSpec() != NULL &&
         *type.getParentTypeSpec() == WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >::getClassTypeSpec() &&
         type.getElementTypeSpec() != NULL &&
         rlType::getClassTypeSpec().isDerivedOf(*type.getElementTypeSpec()) )
    {
        
        // First generate an empty model vector of the desired type
        RevObject *empty_container = Workspace::userWorkspace().makeNewDefaultObject( type.getType() );
        Container *the_converted_container = dynamic_cast<Container*>( empty_container );
        
        for (size_t i=0; i<size(); ++i)
        {
            
            const rlType& orgElement = (*this->value)[ i ];
            if ( type.getElementTypeSpec() != NULL  && orgElement.isType( *type.getElementTypeSpec() ) == true )
            {
                the_converted_container->push_back( orgElement );
            }
//            else
//            {
//                RevObject *convObj = orgElement.convertTo( *type.getElementTypeSpec() );
//                theConvertedContainer->push_back( *convObj );
//                delete convObj;
//            }
        
        }
    
        // return the contain
        return empty_container;
    }
    
    
    // Call the base class if all else fails. This will eventually throw an error if the type conversion is not supported.
    return WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >::convertTo( type );
}



/**
 * Map calls to member methods.
 *
 * @param name method called
 * @param args arguments to the method
 * @param[out] found whether the method has been found already
 */
template <typename rlType>
RevPtr<RevVariable> WorkspaceVector<rlType>::executeMethod( std::string const &name, const std::vector<Argument> &args, bool &found )
{
    
    if ( name == "size" )
    {
        found = true;
        
        // return a new RevVariable with the size of this container
        return RevPtr<RevVariable>( new RevVariable( new Natural( size() ), "" ) );
    }
    else if ( name == "[]" )
    {
        found = true;
        
        long index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() - 1;
        return RevPtr<RevVariable>( new RevVariable( getElement( index ) ) );
    }
    else if ( name == "append" )
    {
        found = true;
        
        const rlType& x = static_cast<const rlType&>( args[0].getVariable()->getRevObject() );
        this->push_back(x);
        
        return NULL;
    }
    else if ( name == "erase" )
    {
        found = true;
        const rlType& x = static_cast<const rlType&>( args[0].getVariable()->getRevObject() );
        RevBayesCore::RbConstIterator<rlType> pos = this->find( x );
        if ( pos != this->end() )
        {
            this->erase( pos );
        }
        
        return NULL;
    }
    
    return WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >::executeMethod( name, args, found );
}



/**
 * Get Rev type of object. This is the Rev object element type followed by
 * a single set of square brackets. This provides a nice and convenient way
 * of specifying generic types of vectors for all Rev object types.
 */
template <typename rlType>
const std::string& WorkspaceVector<rlType>::getClassType(void)
{
    static std::string rev_type = rlType::getClassType() + "[]";
    
    return rev_type;
}


/**
 * Get class type spec describing type of object. Note that we have to
 * use the special version of the TypeSpec constructor because we derive
 * directly from WorkspaceContainer and have the special "rlType[]" type
 * specification.
 */
template <typename rlType>
const TypeSpec& WorkspaceVector<rlType>::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >::getClassTypeSpec(), &rlType::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get element at specified position.
 */
template <typename rlType>
rlType* WorkspaceVector<rlType>::getElement(size_t idx) const
{
    const rlType &element = this->operator[](idx);
    return element.clone();
}


/** Get the type spec of this class. */
template <typename rlType>
const TypeSpec& WorkspaceVector<rlType>::getTypeSpec(void) const
{
    return getClassTypeSpec();  // This should do the trick; there should be a separate version of the function for each template type
}


/**
 * Initialize the methods associated with this object, and their argument rules.
 */
template <typename rlType>
void WorkspaceVector<rlType>::initMethods( void )
{
    
    ArgumentRules* size_arg_rules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), size_arg_rules) );
    
    ArgumentRules* element_arg_rules = new ArgumentRules();
    element_arg_rules->push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The index of the element.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", rlType::getClassTypeSpec(), element_arg_rules ) );

    ArgumentRules* append_arg_rules = new ArgumentRules();
    append_arg_rules->push_back( new ArgumentRule( "x", rlType::getClassTypeSpec(), "The element that you want to add.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "append", RlUtils::Void, append_arg_rules) );

    ArgumentRules* erase_arg_rules = new ArgumentRules();
    erase_arg_rules->push_back( new ArgumentRule( "x", rlType::getClassTypeSpec(), "The element that you want to erase.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "erase", RlUtils::Void, erase_arg_rules) );

}


/**
 * In this function we check whether this type is convertible to some other
 * Rev object type. Here we focus entirely on supporting conversion to
 * other generic vectors with compatible elements. This is not done automatically
 * because of the templating: a vector of RealPos does not inherit from a vector
 * of Real, for example.
 */
template <typename rlType>
double WorkspaceVector<rlType>::isConvertibleTo( const TypeSpec& type, bool once ) const
{

    
    if ( type.getParentTypeSpec() != NULL &&
         *type.getParentTypeSpec() == WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >::getClassTypeSpec() &&
         type.getElementTypeSpec() != NULL &&
         rlType::getClassTypeSpec().isDerivedOf(*type.getElementTypeSpec()) )
    {
        // yes we can
        return 0.0;
    }
    else if ( type == WorkspaceVector<RevObject>::getClassTypeSpec() )
    {
        // yes we can
        return 0.0;
    }
    
    return WorkspaceToCoreWrapperObject<RevBayesCore::RbVector<rlType> >::isConvertibleTo( type, once );
}


/**
 * Drop an element from the back of the vector.
 *
 * @return the dropped element
 */
template <typename rlType>
rlType&  WorkspaceVector<rlType>::pop_back( void )
{
    return this->value->pop_back();
}


/**
 * Drop an element from the front of the vector.
 *
 * @return the dropped element
 */
template <typename rlType>
rlType&  WorkspaceVector<rlType>::pop_front( void )
{
    return this->value->erase( this->value->begin() );
}


/**
 * Push a Rev object element onto the back of the vector.
 *
 * @param x the new element
 */
template <typename rlType>
void WorkspaceVector<rlType>::push_back( const RevObject &x )
{
    // cast the object
    const rlType *x_converted = dynamic_cast< const rlType* >( &x );
    
    if ( x_converted == NULL )
    {
        throw RbException() << "Could not append an element of type " << x.getType() << " to a vector of type " <<  this->getType() ;
    }
    
    // Push it onto the back of the elements vector
    this->value->push_back( *x_converted );
}


/**
 * Push a Rev object element onto the front of the vector.
 *
 * @param x the new element
 */
template <typename rlType>
void WorkspaceVector<rlType>::push_front( const rlType &x )
{
    
    // Push it onto the front of the elements vector
    this->value->insert( this->value->begin(), x );
}


/**
 * Print the value of the vector. Because some workspace objects
 * require a lot of space to print their value, we opt for a format that should
 * work well for such cases.
 */
template<typename rlType>
void WorkspaceVector<rlType>::printValue( std::ostream& o, bool ) const
{
    o << std::endl;

    std::stringstream s;
    if ( size() == 1 )
    {
        s << getClassType() << " vector with 1 value";
    }
    else
    {
        s << getClassType() << " vector with " << size() << " values";
    }
    o << s.str() << std::endl;
    
    for ( size_t i = 0; i < s.str().size(); ++i )
    {
        o << "=";
    }
    o << std::endl << std::endl;

    for ( size_t i = 0; i < this->value->size(); ++i )
    {
        o << "[" << i + 1 << "]" << std::endl;
        RevObject& tmp = (*this->value)[i];
        tmp.printValue(o,true);
        o << std::endl << std::endl;
    }

    o << std::endl;
}



/**
 * Get the size of the vector.
 */
template <typename rlType>
size_t WorkspaceVector<rlType>::size( void ) const
{
    return this->value->size();
}


#endif
