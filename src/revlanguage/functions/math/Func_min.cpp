#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>
#include "MinFunction.h"
#include "Func_min.h"
#include "ModelVector.h"
#include "Real.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_min::Func_min( void ) : TypedFunction<Real>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_min* Func_min::clone( void ) const
{
    
    return new Func_min( *this );
}


RevBayesCore::TypedFunction<double>* Func_min::createFunction( void ) const
{
    
    RevBayesCore::DagNode* node = this->args[0].getVariable()->getRevObject().getDagNode();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* arg_vector = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*>(node);
    RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >* arg_matrix = dynamic_cast<RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >*>(node);

    RevBayesCore::MinFunction* f = arg_vector != NULL ? new RevBayesCore::MinFunction( arg_vector ) : new RevBayesCore::MinFunction( arg_matrix );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_min::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> paramTypes;
        paramTypes.push_back( ModelVector<Real>::getClassTypeSpec() );
        paramTypes.push_back( MatrixReal::getClassTypeSpec() );
        argumentRules.push_back( new ArgumentRule( "x", paramTypes, "A vector of values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_min::getClassType(void)
{
    
    static std::string rev_type = "Func_min";
    
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_min::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_min::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "min";
    
    return f_name;
}


const TypeSpec& Func_min::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
