#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_getContinuousCharacterAsVector.h"
#include "ModelVector.h"
#include "GetContinuousCharacterAsVectorFunction.h"
#include "RlContinuousCharacterData.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "RlString.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class ContinuousCharacterData; }
namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** Default constructor */
Func_getContinuousCharacterAsVector::Func_getContinuousCharacterAsVector( void ) : TypedFunction<ModelVector<Real> >()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_getContinuousCharacterAsVector* Func_getContinuousCharacterAsVector::clone( void ) const
{

    return new Func_getContinuousCharacterAsVector( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_getContinuousCharacterAsVector::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<RevBayesCore::ContinuousCharacterData>* data = static_cast<const ContinuousCharacterData &>( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<std::int64_t>* site_index = static_cast<const Natural &>( args[1].getVariable()->getRevObject() ).getDagNode();

    const std::string& order = static_cast<const RlString &>( args[2].getVariable()->getRevObject() ).getValue();
    RevBayesCore::GetContinuousCharacterAsVectorFunction::VECTOR_ORDER ord;
    if (order == "alphabetical")
    {
        ord = RevBayesCore::GetContinuousCharacterAsVectorFunction::VECTOR_ORDER::ALPHABETICAL;
    }
    else
    {
        throw RbException("argument order_by currently only supports \"alphabetical\"");
    }

    RevBayesCore::GetContinuousCharacterAsVectorFunction* f = new RevBayesCore::GetContinuousCharacterAsVectorFunction( data, site_index, ord );

    return f;
}


/** Get argument rules */
const ArgumentRules& Func_getContinuousCharacterAsVector::getArgumentRules( void ) const
{

    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;

    if ( rules_set == false )
    {

        argument_rules.push_back( new ArgumentRule( "data", ContinuousCharacterData::getClassTypeSpec(), "The character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "variance_site", Natural::getClassTypeSpec(), "The site that stores the trait variance.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        std::vector<std::string> vectorOrderTypes;
        vectorOrderTypes.push_back( "alphabetical" );
        argument_rules.push_back( new OptionRule ("order_by", new RlString("alphabetical"), vectorOrderTypes, "The order of the elements in the vector. Option \"alphabetical\" returns the the character ordered alphbetically by the species name.") );

        rules_set = true;
    }

    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_getContinuousCharacterAsVector::getClassType(void)
{

    static std::string rev_type = "Func_getContinuousCharacterAsVector";

    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_getContinuousCharacterAsVector::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_getContinuousCharacterAsVector::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnGetContinuousCharacterAsVector";

    return f_name;
}

std::vector<std::string> Func_getContinuousCharacterAsVector::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "fnGetContChar" );

    return a_names;
}

/** Get type spec */
const TypeSpec& Func_getContinuousCharacterAsVector::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}
