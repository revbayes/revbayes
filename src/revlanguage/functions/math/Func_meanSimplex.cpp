#include <iosfwd>
#include <string>
#include <vector>

#include "MeanFunction.h"
#include "Func_meanSimplex.h"
#include "RlSimplex.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

/** default constructor */
Func_meanSimplex::Func_meanSimplex( void ) : TypedFunction<RealPos>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_meanSimplex* Func_meanSimplex::clone( void ) const
{
    
    return new Func_meanSimplex( *this );
}


RevBayesCore::TypedFunction<double>* Func_meanSimplex::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* arg = static_cast<RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* >( this->args[0].getVariable()->getRevObject().getDagNode() );
    RevBayesCore::MeanFunction* f = new RevBayesCore::MeanFunction( arg );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_meanSimplex::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x", Simplex::getClassTypeSpec(), "A vector of numbers.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_meanSimplex::getClassType(void)
{
    
    static std::string rev_type = "Func_meanSimplex";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_meanSimplex::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_meanSimplex::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "mean";
    
    return f_name;
}


const TypeSpec& Func_meanSimplex::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
