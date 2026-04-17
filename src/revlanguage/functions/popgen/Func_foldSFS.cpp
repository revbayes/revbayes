#include "Func_foldSFS.h"

#include <string>
#include <vector>

#include "FoldSFSFunction.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "Natural.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_foldSFS::Func_foldSFS(void) : TypedFunction< ModelVector< Natural > >()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objects.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the function.
 */
Func_foldSFS* Func_foldSFS::clone( void ) const
{
    return new Func_foldSFS( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RbVector<std::int64_t> >* Func_foldSFS::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode< RevBayesCore::RbVector<std::int64_t> >* sfs =
        static_cast<const ModelVector<Natural> &>( args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::FoldSFSFunction* f = new RevBayesCore::FoldSFSFunction( sfs );

    return f;
}


const ArgumentRules& Func_foldSFS::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "sfs", ModelVector<Natural>::getClassTypeSpec(), "An unfolded site frequency spectrum (length n+1 for n individuals).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_foldSFS::getClassType(void)
{
    static std::string rev_type = "Func_foldSFS";
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_foldSFS::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_foldSFS::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnFoldSFS";
    
    return f_name;
}


const TypeSpec& Func_foldSFS::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
