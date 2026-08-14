#include "Func_getNumProcesses.h"

#include "Natural.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevLanguage;

/** Default constructor */
Func_getNumProcesses::Func_getNumProcesses( void ) : Procedure()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objects.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the function.
 */
Func_getNumProcesses* Func_getNumProcesses::clone( void ) const
{

    return new Func_getNumProcesses( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_getNumProcesses::execute( void )
{

    int num_processes = 1;

#ifdef RB_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
#endif

    return new RevVariable( new Natural( num_processes ) );
}


/** Get argument rules */
const ArgumentRules& Func_getNumProcesses::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();

    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_getNumProcesses::getClassType(void)
{

    static std::string rev_type = "Func_getNumProcesses";

    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_getNumProcesses::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_getNumProcesses::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "getNumProcesses";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_getNumProcesses::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get return type */
const TypeSpec& Func_getNumProcesses::getReturnType( void ) const
{

    static TypeSpec return_typeSpec = Natural::getClassTypeSpec();

    return return_typeSpec;
}
