#include <sstream>
#include <string>

#include "Natural.h"
#include "Probability.h"
#include "RealPos.h"
#include "RbException.h"
#include "TypeSpec.h"
#include "Real.h"
#include "RevObject.h"
#include "TypedDagNode.h"

using namespace RevLanguage;
    
/** Default constructor */
RealPos::RealPos( void ) : Real( 1.0 )
{

}


/** Construct from double */
RealPos::RealPos( RevBayesCore::TypedDagNode<double> *x ) : Real( x )
{
    
    if ( x->getValue() < 0.0 )
    {
        throw RbException( "Negative value for " + getClassType() );
    }
    
}


/** Construct from double */
RealPos::RealPos( double x ) : Real( x )
{

    if ( x < 0.0 )
    {
        throw RbException( "Negative value for " + getClassType() );
    }
    
}


/** Construct from int */
RealPos::RealPos( long x ) : Real( double(x) )
{

    if ( x < 0 )
    {
        throw RbException( "Negative value for " + getClassType() );
    }
    
}

/** Construct from bool */
RealPos::RealPos( bool x ) : Real( double(x) ) {}


/**
 * Generic addition operator.
 * We test if the rhs is of a type that we use for a specialized addition operation.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RevObject* RealPos::add( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
    {
        return add( static_cast<const RealPos&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf( Natural::getClassTypeSpec() ) )
    {
        return add( static_cast<const Natural&>( rhs ) );
    }
    
    return Real::add( rhs );
}


/**
 * Specialized addition operation between two positive real numbers.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RealPos* RealPos::add(const RevLanguage::Natural &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() + rhs.getValue() );
    
    return n;
}


/**
 * Specialized addition operation between two positive real numbers.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RealPos* RealPos::add(const RevLanguage::RealPos &rhs) const
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
RealPos* RealPos::clone( void ) const
{

	return new RealPos( *this );
}

RevObject* RealPos::convertTo( const TypeSpec& type ) const
{
    
    if ( type == Real::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<RealPos,Real>(this);
    if ( type == Probability::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<RealPos,Probability>(this);

    return Real::convertTo( type );
}


/**
 * Generic division operator.
 * We test if the rhs is of a type that we use for a specialized division operation.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RevObject* RealPos::divide( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
    {
        return divide( static_cast<const RealPos&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf( Natural::getClassTypeSpec() ) )
    {
        return divide( static_cast<const Natural&>( rhs ) );
    }
    
    return Real::divide( rhs );
}


/**
 * Specialized division operation between two positive real numbers.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RealPos* RealPos::divide(const RevLanguage::Natural &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() / rhs.getValue() );
    
    return n;
}


/**
 * Specialized division operation between two positive real numbers.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RealPos* RealPos::divide(const RevLanguage::RealPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() / rhs.getValue() );
    
    return n;
}


/** Get Rev type of object */
const std::string& RealPos::getClassType(void)
{
    
    static std::string rev_type = "RealPos";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& RealPos::getClassTypeSpec(void) { 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Real::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& RealPos::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Is convertible to type? */
double RealPos::isConvertibleTo(const TypeSpec& type, bool once) const
{
    
    if ( type == Real::getClassTypeSpec() )
    {
        return 0.2;
    }
    else if ( once == true && type == Probability::getClassTypeSpec() && dag_node->getValue() <= 1.0 )
    {
        return 0.1;
    }
    
    double tmp = Real::isConvertibleTo(type, once);
    return ( (tmp == -1.0) ? -1.0 : (tmp+0.2));
}


/**
 * Generic multiplication operator.
 * We test if the rhs is of a type that we use for a specialized multiplication operation.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
RevObject* RealPos::multiply( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( RealPos::getClassTypeSpec() ) )
    {
        return multiply( static_cast<const RealPos&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf( Natural::getClassTypeSpec() ) )
    {
        return multiply( static_cast<const Natural&>( rhs ) );
    }
    
    return Real::multiply( rhs );
}


/**
 * Specialized multiplication operation between two positive real numbers.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
RealPos* RealPos::multiply(const RevLanguage::Natural &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() * rhs.getValue() );
    
    return n;
}


/**
 * Specialized multiplication operation between two positive real numbers.
 * The return value is also of type positive real number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
RealPos* RealPos::multiply(const RevLanguage::RealPos &rhs) const
{
    
    RealPos *n = new RealPos( dag_node->getValue() * rhs.getValue() );
    
    return n;
}
