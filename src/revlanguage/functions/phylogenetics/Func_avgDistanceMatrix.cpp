#include "Func_avgDistanceMatrix.h"

#include "AverageDistanceMatrix.h"
#include "DistanceMatrix.h"
#include "ModelVector.h"
#include "RlDeterministicNode.h"
#include "AvgDistanceMatrixFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "RlMatrixBoolean.h"
#include "RlDistanceMatrix.h"
#include "TypeSpec.h"

namespace RevBayesCore { class DistanceMatrix; }

using namespace RevLanguage;


/** default constructor */
Func_avgDistanceMatrix::Func_avgDistanceMatrix( void ) : TypedFunction< AverageDistanceMatrix >( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_avgDistanceMatrix* Func_avgDistanceMatrix::clone( void ) const
{
    
    return new Func_avgDistanceMatrix( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::AverageDistanceMatrix>* Func_avgDistanceMatrix::createFunction( void ) const
{
    
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::DistanceMatrix> >* matrixVector = static_cast<const ModelVector<DistanceMatrix> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    RevBayesCore::AvgDistanceMatrixFunction * f = new RevBayesCore::AvgDistanceMatrixFunction( matrixVector );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_avgDistanceMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set     = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "distanceMatrices", ModelVector<DistanceMatrix>::getClassTypeSpec(), "The vector of distance matrices.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_avgDistanceMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_avgDistanceMatrix";
    
    return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_avgDistanceMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedFunction<AverageDistanceMatrix>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_avgDistanceMatrix::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "avgDistanceMatrix" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_avgDistanceMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instances of this class
    std::string f_name = "fnAverageDistanceMatrix";
    
    return f_name;
}


const TypeSpec& Func_avgDistanceMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

