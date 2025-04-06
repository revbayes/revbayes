#include <sstream>
#include <string>

#include "RlBoolean.h"
#include "Natural.h"
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
Natural::Natural( void ) : Integer( 0 )
{

}


Natural::Natural( RevBayesCore::TypedDagNode<std::int64_t> *v ) : Integer( v )
{
    
}


/** Construct from int */
Natural::Natural( std::int64_t x ) : Integer( x )
{

    if ( x < 0 )
    {
        throw RbException() << "Negative value for " << getClassType() ; 
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
RevObject* Natural::add( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( Natural::getClassTypeSpec() ) )
    {
        return add( static_cast<const Natural&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
    {
        return add( static_cast<const RealPos&>( rhs ) );
    }
    
    return Integer::add( rhs );
}


/**
 * Specialized addition operation between two natural numbers.
 * The return value is also of type natural number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
Natural* Natural::add(const RevLanguage::Natural &rhs) const
{
    
    Natural *n = new Natural( dag_node->getValue() + rhs.getValue() );
    
    return n;
}


/**
 * Specialized addition operation between a natural and a positiove real number.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RealPos* Natural::add(const RevLanguage::RealPos &rhs) const
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
Natural* Natural::clone( void ) const
{

	return new Natural( *this );
}


/** Convert to type. The caller manages the returned object. */
RevObject* Natural::convertTo( const TypeSpec& type ) const
{
    if ( type == RlBoolean::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<Natural,RlBoolean>(this);
    
    if ( type == Real::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<Natural,Real>(this);   
    if ( type == RealPos::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<Natural,RealPos>(this);
    
    if ( type == Probability::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<Natural,Probability>(this);

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
RevObject* Natural::divide( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( Natural::getClassTypeSpec() ) )
    {
        return divide( static_cast<const Natural&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
    {
        return divide( static_cast<const RealPos&>( rhs ) );
    }
    
    return Integer::divide( rhs );
}


/**
 * Specialized division operation between two natural numbers.
 * The return value is of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RealPos* Natural::divide(const RevLanguage::Natural &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() / double( rhs.getValue() ) );
    
    return n;
}


/**
 * Specialized division operation between a natural and a positiove real number.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RealPos* Natural::divide(const RevLanguage::RealPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() / rhs.getValue() );
    
    return n;
}


/** Get Rev type of object */
const std::string& Natural::getClassType(void)
{
    
    static std::string rev_type = "Natural";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Natural::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Integer::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& Natural::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Is convertible to type? */
double Natural::isConvertibleTo( const TypeSpec& type, bool once ) const
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
RevObject* Natural::multiply( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( Natural::getClassTypeSpec() ) )
        return multiply( static_cast<const Natural&>( rhs ) );
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
        return multiply( static_cast<const RealPos&>( rhs ) );
    
    return Integer::multiply( rhs );
}


/**
 * Specialized multiplication operation between two natural numbers.
 * The return value is also of type natural number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
Natural* Natural::multiply(const RevLanguage::Natural &rhs) const
{
    
    Natural *n = new Natural( dag_node->getValue() * rhs.getValue() );
    
    return n;
}


/**
 * Specialized multiplication operation between a natural and a positive real number.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
RealPos* Natural::multiply(const RevLanguage::RealPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() * rhs.getValue() );
    
    return n;
}
