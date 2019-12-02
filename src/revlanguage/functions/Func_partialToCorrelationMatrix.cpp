//
//  Func_varianceCovarianceMatrix.cpp
//  RevBayesCore
//
//  Created by Michael R. May on 08/29/17.
//

#include "Func_partialToCorrelationMatrix.h"
#include "MatrixReal.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "PartialToCorrelationFunction.h"

using namespace RevLanguage;


/** default constructor */
Func_partialToCorrelationMatrix::Func_partialToCorrelationMatrix(void) : TypedFunction<MatrixRealSymmetric>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_partialToCorrelationMatrix* Func_partialToCorrelationMatrix::clone(void) const {
    
    return new Func_partialToCorrelationMatrix(*this);
}


RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_partialToCorrelationMatrix::createFunction(void) const {

    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* r = static_cast<const MatrixReal& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::PartialToCorrelation* f = new RevBayesCore::PartialToCorrelation( r );
    return f;
    
}


/* Get argument rules */
const ArgumentRules& Func_partialToCorrelationMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "partialCorrelations", MatrixReal::getClassTypeSpec(), "A matrix of partial correlation coefficients.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    return argumentRules;
}


const std::string& Func_partialToCorrelationMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_partialToCorrelationMatrix";
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_partialToCorrelationMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_partialToCorrelationMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPartialToCorr";
    
    return f_name;
}


const TypeSpec& Func_partialToCorrelationMatrix::getTypeSpec(void) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
