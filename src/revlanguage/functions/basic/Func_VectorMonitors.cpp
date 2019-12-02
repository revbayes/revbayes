#include "ArgumentRule.h"
#include "Func_VectorMonitors.h"
#include "RlMonitor.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"


using namespace RevLanguage;


/** Default constructor */
Func_VectorMonitors::Func_VectorMonitors( void ) : Procedure()
{
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_VectorMonitors* Func_VectorMonitors::clone( void ) const
{
    return new Func_VectorMonitors( *this );
}


/** Execute function: Construct an empty VectorMonitors. */
RevPtr<RevVariable> Func_VectorMonitors::execute( void )
{
    
    WorkspaceVector<Monitor> *vec = new WorkspaceVector<Monitor>;
    RevVariable *var = new RevVariable(vec);
    
    return RevPtr<RevVariable>( var );
}


/** Get argument rules */
const ArgumentRules& Func_VectorMonitors::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object (static) */
const std::string& Func_VectorMonitors::getClassType( void )
{
    static std::string rev_type = "Func_VectorMonitors";
    
    return rev_type;
}


/** Get Rev type spec of object (static) */
const TypeSpec& Func_VectorMonitors::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_VectorMonitors::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "VectorMonitors";
    
    return f_name;
}

const TypeSpec& Func_VectorMonitors::getReturnType(void) const
{
    static TypeSpec return_typeSpec = WorkspaceVector<Monitor>::getClassTypeSpec();
    return return_typeSpec;
}


/** Get Rev type spec of object from an instance */
const TypeSpec& Func_VectorMonitors::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

