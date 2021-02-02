#include <iosfwd>

#include "Func_VectorMoves.h"
#include "RlMove.h"
#include "TypeSpec.h"
#include "WorkspaceVector.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "WorkspaceToCoreWrapperObject.h"


using namespace RevLanguage;


/** Default constructor */
Func_VectorMoves::Func_VectorMoves( void ) : Procedure()
{
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_VectorMoves* Func_VectorMoves::clone( void ) const
{
    return new Func_VectorMoves( *this );
}


/** Execute function: Construct an empty VectorMoves. */
RevPtr<RevVariable> Func_VectorMoves::execute( void )
{
    
    WorkspaceVector<Move> *vec = new WorkspaceVector<Move>;
    RevVariable *var = new RevVariable(vec);
    
    return RevPtr<RevVariable>( var );
}


/** Get argument rules */
const ArgumentRules& Func_VectorMoves::getArgumentRules( void ) const
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
const std::string& Func_VectorMoves::getClassType( void )
{
    static std::string rev_type = "Func_VectorMoves";
    
    return rev_type;
}


/** Get Rev type spec of object (static) */
const TypeSpec& Func_VectorMoves::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_VectorMoves::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "VectorMoves";
    
    return f_name;
}

const TypeSpec& Func_VectorMoves::getReturnType(void) const
{
    static TypeSpec return_typeSpec = WorkspaceVector<Move>::getClassTypeSpec();
    return return_typeSpec;
}


/** Get Rev type spec of object from an instance */
const TypeSpec& Func_VectorMoves::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

