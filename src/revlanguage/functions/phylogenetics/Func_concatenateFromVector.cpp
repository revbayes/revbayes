#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_concatenateFromVector.h"
#include "OptionRule.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Procedure.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlFunction.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** Default constructor */
Func_concatenateFromVector::Func_concatenateFromVector( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_concatenateFromVector* Func_concatenateFromVector::clone( void ) const
{
    
    return new Func_concatenateFromVector( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_concatenateFromVector::execute( void )
{
    const ModelVector<AbstractHomologousDiscreteCharacterData>& v = static_cast<const ModelVector<AbstractHomologousDiscreteCharacterData> &>( args[0].getVariable()->getRevObject() );
    const std::string& type = static_cast<const RlString &>( args[1].getVariable()->getRevObject() ).getValue();
    
    AbstractHomologousDiscreteCharacterData* a = dynamic_cast<const AbstractHomologousDiscreteCharacterData *>( v.getElement(0) )->clone();
    for (size_t i = 1; i < v.size(); ++i)
    {
        const AbstractHomologousDiscreteCharacterData* b = dynamic_cast<const AbstractHomologousDiscreteCharacterData *>( v.getElement(i) );
        
        a->concatenate( *b, type );
        
    }
    
    return new RevVariable( a );
}


/** Get argument rules */
const ArgumentRules& Func_concatenateFromVector::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "x", ModelVector<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec(), "First character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<std::string> optionsCondition;
        optionsCondition.push_back( "union" );
        optionsCondition.push_back( "intersection" );
        argument_rules.push_back( new OptionRule( "merge", new RlString("Taxa must match"), optionsCondition, "How to merge differing taxa" ) );
        
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_concatenateFromVector::getClassType(void)
{
    
    static std::string rev_type = "Func_concatenateFromVector";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_concatenateFromVector::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_concatenateFromVector::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "concatenate";
    
    return f_name;
}

/** Get type spec */
const TypeSpec& Func_concatenateFromVector::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_concatenateFromVector::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlBoolean::getClassTypeSpec();
    
    return return_typeSpec;
}

