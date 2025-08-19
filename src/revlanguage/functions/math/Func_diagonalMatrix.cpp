#include "Func_diagonalMatrix.h"

#include <cstddef>

#include "Natural.h"
#include "RlMatrixRealSymmetric.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MatrixReal.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** default constructor */
Func_diagonalMatrix::Func_diagonalMatrix( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_diagonalMatrix* Func_diagonalMatrix::clone( void ) const
{
    
    return new Func_diagonalMatrix( *this );
}


RevPtr<RevVariable> Func_diagonalMatrix::execute( void )
{
    
    int n = (int)static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
    RevBayesCore::MatrixReal *m = new RevBayesCore::MatrixReal(n, n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        (*m)[i][i] = 1.0;
    }
    
    return new RevVariable( new MatrixRealSymmetric( m ) );
}



/* Get argument rules */
const ArgumentRules& Func_diagonalMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule("n", Natural::getClassTypeSpec(), "The number of rows/columns (dimension).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_diagonalMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_diag";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_diagonalMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_diagonalMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "diagonalMatrix";
    
    return f_name;
}


const TypeSpec& Func_diagonalMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_diagonalMatrix::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = MatrixRealSymmetric::getClassTypeSpec();
    return return_typeSpec;
}
