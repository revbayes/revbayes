#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "Integer.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "Real.h"
#include "RealPos.h"
#include "TypeSpec.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
RlBoolean::RlBoolean(void) : ModelObject<RevBayesCore::Boolean>()
{

}

/** Construct from bool */
RlBoolean::RlBoolean(RevBayesCore::TypedDagNode<RevBayesCore::Boolean> *v) : ModelObject<RevBayesCore::Boolean>( v )
{
    
}

/** Construct from bool */
RlBoolean::RlBoolean(bool v) : ModelObject<RevBayesCore::Boolean>( new RevBayesCore::Boolean(v) )
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
RlBoolean* RlBoolean::clone(void) const {

	return new RlBoolean(*this);
}


/** Convert to type. The caller manages the returned object. */
RevObject* RlBoolean::convertTo(const TypeSpec& type) const
{

    if (type == Natural::getClassTypeSpec()) return RlUtils::RlTypeConverter::convertTo<RlBoolean,Natural>(this);  
    if (type == Integer::getClassTypeSpec()) return RlUtils::RlTypeConverter::convertTo<RlBoolean,Integer>(this);

    if (type == Real::getClassTypeSpec()) return RlUtils::RlTypeConverter::convertTo<RlBoolean,Real>(this);
    if (type == RealPos::getClassTypeSpec()) return RlUtils::RlTypeConverter::convertTo<RlBoolean,RealPos>(this);

    return RevObject::convertTo(type);
}


/** Get Rev type of object */
const std::string& RlBoolean::getClassType(void)
{
    
    static std::string rev_type = "Bool";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& RlBoolean::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& RlBoolean::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Is convertible to type? */
double RlBoolean::isConvertibleTo(const TypeSpec& type, bool once) const
{

    if ( type == Natural::getClassTypeSpec() )
    {
        return 0.1;
    }
    else if ( type == Integer::getClassTypeSpec() )
    {
        return 0.3;
    }
    else if ( type == RealPos::getClassTypeSpec() )
    {
        return 0.2;
    }
    else if ( type == Real::getClassTypeSpec() )
    {
        return 0.4;
    }
    
    return RevObject::isConvertibleTo(type, once);
}


/** Print value for user */
void RlBoolean::printValue(std::ostream &o) const
{

    o << (dag_node->getValue() ? "true" : "false");
}

