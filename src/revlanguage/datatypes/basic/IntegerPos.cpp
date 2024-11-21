#include <sstream>
#include <cstdint>
#include <string>

#include "RlBoolean.h"
#include "IntegerPos.h"
#include "Integer.h"
#include "Probability.h"
#include "RealPos.h"
#include "Real.h"
#include "RbException.h"
#include "RlDiscreteCharacterState.h"
#include "RlString.h"
#include "StandardState.h"
#include "TypeSpec.h"
#include "RevObject.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** Default constructor */
IntegerPos::IntegerPos( void ) : Natural( 0L )
{

}


IntegerPos::IntegerPos( RevBayesCore::TypedDagNode<std::int64_t> *v ) : Natural( v )
{
    
}


/** Construct from Natural */
IntegerPos::IntegerPos( std::int64_t x ) : Natural( x )
{

    if ( x < 1 )
    {
        throw RbException( "Negative or zero value for " + getClassType() );
    }
    
}


/**
 * Generic addition operator.
 * We test if the rhs is of a type that we use for a specialized addition operation.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RevObject* IntegerPos::add( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( IntegerPos::getClassTypeSpec() ) )
    {
        return add( static_cast<const IntegerPos&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
    {
        return add( static_cast<const RealPos&>( rhs ) );
    }
    
    return Integer::add( rhs );
}


/**
 * Specialized addition operation between two IntegerPos numbers.
 * The return value is also of type IntegerPos number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
IntegerPos* IntegerPos::add(const RevLanguage::IntegerPos &rhs) const
{
    
    IntegerPos *n = new IntegerPos( dag_node->getValue() + rhs.getValue() );
    
    return n;
}


/**
 * Specialized addition operation between a positive integer and a positive real number.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RealPos* IntegerPos::add(const RevLanguage::RealPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() + rhs.getValue() );
    
    return n;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
IntegerPos* IntegerPos::clone( void ) const
{

	return new IntegerPos( *this );
}


/** Convert to type. The caller manages the returned object. */
RevObject* IntegerPos::convertTo( const TypeSpec& type ) const
{
    if ( type == RlBoolean::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<IntegerPos,RlBoolean>(this);
    
    if ( type == Real::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<IntegerPos,Real>(this);   
    if ( type == RealPos::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<IntegerPos,RealPos>(this);
    
    if ( type == Probability::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<IntegerPos,Probability>(this);

    if ( type == RlString::getClassTypeSpec() )
    {
        std::ostringstream o;
        printValue( o, true );
        return new RlString( o.str() );
    }
    
    if ( type == DiscreteCharacterState::getClassTypeSpec() )
    {        
        std::ostringstream o;
        printValue( o, true );
        return new DiscreteCharacterState( RevBayesCore::StandardState( o.str() ) );
    }

    return Integer::convertTo( type );
}


/**
 * Generic division operator.
 * We test if the rhs is of a type that we use for a specialized division operation.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RevObject* IntegerPos::divide( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( IntegerPos::getClassTypeSpec() ) )
    {
        return divide( static_cast<const IntegerPos&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
    {
        return divide( static_cast<const RealPos&>( rhs ) );
    }
    
    return IntegerPos::divide( rhs );
}


/**
 * Specialized division operation between two positive integer numbers.
 * The return value is of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RealPos* IntegerPos::divide(const RevLanguage::IntegerPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() / double( rhs.getValue() ) );
    
    return n;
}


/**
 * Specialized division operation between a positive integer and a positive real number.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RealPos* IntegerPos::divide(const RevLanguage::RealPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() / rhs.getValue() );
    
    return n;
}


/** Get Rev type of object */
const std::string& IntegerPos::getClassType(void)
{
    
    static std::string rev_type = "IntegerPos";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& IntegerPos::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Integer::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& IntegerPos::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Is convertible to type? */
double IntegerPos::isConvertibleTo( const TypeSpec& type, bool once ) const
{

    if ( type == RlBoolean::getClassTypeSpec() )
    {
        return 0.5;
    }
    
    if ( type == Real::getClassTypeSpec() )
    {
        return 0.3;
    }
    
    if ( type == RealPos::getClassTypeSpec() )
    {
        return 0.2;
    }
    
    if ( once == true && type == Probability::getClassTypeSpec() && dag_node->getValue() <= 1 )
    {
        return 0.1;
    }
    
    if ( type == RlString::getClassTypeSpec() )
    {
        return 0.6;
    }
    
    if ( type == DiscreteCharacterState::getClassTypeSpec() )
    {
        return 0.7;
    }
    
    return Integer::isConvertibleTo( type, once );
}


/**
 * Generic multiplication operator.
 * We test if the rhs is of a type that we use for a specialized multiplication operation.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
RevObject* IntegerPos::multiply( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( IntegerPos::getClassTypeSpec() ) )
        return multiply( static_cast<const IntegerPos&>( rhs ) );
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
        return multiply( static_cast<const RealPos&>( rhs ) );
    
    return Integer::multiply( rhs );
}


/**
 * Specialized multiplication operation between two positive integer numbers.
 * The return value is also of type positive integer number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
IntegerPos* IntegerPos::multiply(const RevLanguage::IntegerPos &rhs) const
{
    
    IntegerPos *n = new IntegerPos( dag_node->getValue() * rhs.getValue() );
    
    return n;
}


/**
 * Specialized multiplication operation between a positive integer and a positive real number.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
RealPos* IntegerPos::multiply(const RevLanguage::RealPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() * rhs.getValue() );
    
    return n;
}
