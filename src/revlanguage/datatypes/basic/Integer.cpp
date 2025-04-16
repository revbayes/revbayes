#include <sstream>
#include <string>
#include <vector>

#include "IntegerPos.h"
#include "RlBoolean.h"
#include "Integer.h"
#include "Natural.h"
#include "Probability.h"
#include "RealPos.h"
#include "Real.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/* Default constructor */
Integer::Integer(void) : ModelObject<std::int64_t>()
{
    
}


Integer::Integer( RevBayesCore::TypedDagNode<std::int64_t> *v ) : ModelObject<std::int64_t>( v )
{
    
}



/* Construct from int */
Integer::Integer(std::int64_t v) : ModelObject<std::int64_t>( new std::int64_t(v) )
{

}

/**
 * Generic addition operator.
 * We test if the rhs is of a type that we use for a specialized addition operation.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
RevObject* Integer::add( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( Real::getClassTypeSpec() ) )
        return add( static_cast<const Real&>( rhs ) );
    
    if ( rhs.getTypeSpec().isDerivedOf(  Integer::getClassTypeSpec() ) )
        return add( static_cast<const Integer&>( rhs ) );
    
    return ModelObject<std::int64_t>::add( rhs );
}


/**
 * Specialized addition operation between two real numbers.
 * The return value is also of type real number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
Real* Integer::add(const Real &rhs) const
{
    
    Real *n = new Real( dag_node->getValue() + rhs.getValue() );
    
    return n;
}


/**
 * Specialized addition operation between two real numbers.
 * The return value is also of type real number.
 *
 * \param[in]   rhs     The right hand side operand of the addition operation.
 *
 * \return              A new object holding the sum.
 */
Integer* Integer::add(const Integer &rhs) const
{
    
    Integer *n = new Integer( dag_node->getValue() + rhs.getValue() );
    
    return n;
}


/** 
 * Clone object 
 */
Integer* RevLanguage::Integer::clone(void) const
{

	return new Integer(*this);
}


/** 
 * Convert to type. The caller manages the returned object.
 * Calls the templated convertTo<> function based on the type to return.
 */
RevObject* Integer::convertTo( const TypeSpec& type ) const
{   
    if ( type == RlBoolean::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<Integer,RlBoolean>(this);
    if ( type == Real::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<Integer,Real>(this);

    if ( type == RlString::getClassTypeSpec() ) {
        std::ostringstream o;
        printValue( o, true );
        return new RlString( o.str() );
    }

    if ( type == RealPos::getClassTypeSpec() && dag_node->getValue() >= 0 ) return RlUtils::RlTypeConverter::convertTo<Integer,RealPos>(this);
    if ( type == IntegerPos::getClassTypeSpec() && dag_node->getValue() > 0) return RlUtils::RlTypeConverter::convertTo<Integer,IntegerPos>(this);
    if ( type == Natural::getClassTypeSpec() && dag_node->getValue() >= 0) return RlUtils::RlTypeConverter::convertTo<Integer,Natural>(this);
    
    if ( type == Probability::getClassTypeSpec() ) return RlUtils::RlTypeConverter::convertTo<Integer,Probability>(this);
    
    return RevObject::convertTo( type );
}

/**
  * Specialized decrement operation.
  */
void Integer::decrement( void ) 
{
    
    dag_node->getValue()--;

}


/**
 * Generic division operator.
 * We test if the rhs is of a type that we use for a specialized division operation.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
RevObject* Integer::divide( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( Real::getClassTypeSpec() ) )
        return divide( static_cast<const Real&>( rhs ) );
    
    if ( rhs.getTypeSpec().isDerivedOf(  Integer::getClassTypeSpec() ) )
        return divide( static_cast<const Integer&>( rhs ) );
    
    return ModelObject<std::int64_t>::divide( rhs );
}


/**
 * Specialized division operation between two real numbers.
 * The return value is also of type real number.
 *
 * \param[in]   rhs     The right hand side operand of the divsion operation.
 *
 * \return              A new object holding the ratio.
 */
Real* Integer::divide(const Real &rhs) const
{
    
    Real *n = new Real( dag_node->getValue() / rhs.getValue() );
    
    return n;
}


/**
 * Specialized division operation between two integer numbers.
 * The return value is also of type integer number.
 *
 * \param[in]   rhs     The right hand side operand of the division operation.
 *
 * \return              A new object holding the ratio.
 */
Real* Integer::divide(const Integer &rhs) const
{
    
    Real *n = new Real( dag_node->getValue() / double( rhs.getValue() ) );
    
    return n;
}


/** 
 * Get Rev type of object 
 */
const std::string& Integer::getClassType(void)
{
    
    static std::string rev_type = "Integer";
    
	return rev_type; 
}

/** 
 * Get class type spec describing type of object 
 */
const TypeSpec& Integer::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** 
 * Get type spec 
 */
const TypeSpec& Integer::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Specialized increment operation.
 */
void Integer::increment( void ) 
{
    
    dag_node->getValue()++;
    
}



/** 
 * Is convertible to language object of type? 
 */
double Integer::isConvertibleTo( const TypeSpec& type, bool convert_by_value ) const
{

    if ( type == RlBoolean::getClassTypeSpec() )
    {
        return 0.6;
    }
    
    if ( type == Real::getClassTypeSpec() )
    {
        return 0.4;
    }
    
    if ( type == RlString::getClassTypeSpec() )
    {
        return 0.5;
    }
    
    if ( convert_by_value && type == RealPos::getClassTypeSpec() && dag_node->getValue() >= 0 )
    {
        return 0.3;
    }
    if ( convert_by_value && type == IntegerPos::getClassTypeSpec() && dag_node->getValue() > 0 )
    {
        return 0.1;
    }

    if ( convert_by_value && type == Natural::getClassTypeSpec() && dag_node->getValue() >= 0 )
    {
        return 0.1;
    }
    
    if ( convert_by_value == true && type == Probability::getClassTypeSpec() && dag_node->getValue() <= 1 && dag_node->getValue() >= 0)
    {
        return 0.2;
    }
    
    return RevObject::isConvertibleTo( type, convert_by_value );
}


/**
 * Generic multiplication operator.
 * We test if the rhs is of a type that we use for a specialized multiplication operation.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
RevObject* Integer::multiply( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( Real::getClassTypeSpec() ) )
    {
        return multiply( static_cast<const Real&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf(  Integer::getClassTypeSpec() ) )
    {
        return multiply( static_cast<const Integer&>( rhs ) );
    }
    
    return ModelObject<std::int64_t>::multiply( rhs );
}


/**
 * Specialized multiplication operation between an integer and a real number.
 * The return value is also of type real number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
Real* Integer::multiply(const Real &rhs) const
{
    
    Real *n = new Real( dag_node->getValue() * rhs.getValue() );
    
    return n;
}


/**
 * Specialized multiplication operation between two integer numbers.
 * The return value is also of type integer number.
 *
 * \param[in]   rhs     The right hand side operand of the multiplication operation.
 *
 * \return              A new object holding the product.
 */
Integer* Integer::multiply(const Integer &rhs) const
{
    
    Integer *n = new Integer( dag_node->getValue() * rhs.getValue() );
    
    return n;
}


/**
 * Generic subtraction operator.
 * We test if the rhs is of a type that we use for a specialized subtraction operation.
 *
 * \param[in]   rhs     The right hand side operand of the subtraction operation.
 *
 * \return              A new object holding the difference.
 */
RevObject* Integer::subtract( const RevObject& rhs ) const 
{
    
    if ( rhs.getTypeSpec().isDerivedOf( Real::getClassTypeSpec() ) )
    {
        return subtract( static_cast<const Real&>( rhs ) );
    }
    
    if ( rhs.getTypeSpec().isDerivedOf(  Integer::getClassTypeSpec() ) )
    {
        return subtract( static_cast<const Integer&>( rhs ) );
    }
    
    return ModelObject<std::int64_t>::subtract( rhs );
}


/**
 * Specialized subtraction operation between an integer and a real number.
 * The return value is also of type real number.
 *
 * \param[in]   rhs     The right hand side operand of the subtraction operation.
 *
 * \return              A new object holding the difference.
 */
Real* Integer::subtract(const Real &rhs) const
{
    
    Real *n = new Real( dag_node->getValue() - rhs.getValue() );
    
    return n;
}


/**
 * Specialized subtraction operation between two integer numbers.
 * The return value is also of type ineteger number.
 *
 * \param[in]   rhs     The right hand side operand of the subtraction operation.
 *
 * \return              A new object holding the difference.
 */
Integer* Integer::subtract(const Integer &rhs) const
{
    
    Integer *n = new Integer( dag_node->getValue() - rhs.getValue() );
    
    return n;
}
