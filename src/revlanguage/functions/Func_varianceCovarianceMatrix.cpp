//
//  Func_varianceCovarianceMatrix.cpp
//  RevBayesCore
//
//  Created by John Huelsenbeck on 11/17/14.
//

#include "Func_varianceCovarianceMatrix.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "VarianceCovarianceFunction.h"

using namespace RevLanguage;


/** default constructor */
Func_varianceCovarianceMatrix::Func_varianceCovarianceMatrix(void) : TypedFunction<MatrixRealSymmetric>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_varianceCovarianceMatrix* Func_varianceCovarianceMatrix::clone(void) const
{
    
    return new Func_varianceCovarianceMatrix(*this);
}


RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_varianceCovarianceMatrix::createFunction(void) const {

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sd = static_cast<const ModelVector<RealPos>& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* ce = static_cast<const ModelVector<Real>& >( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    double n_sd = sd->getValue().size();
    double n_ce = ce->getValue().size();
    if ( n_ce != (n_sd * (n_sd-1.0) / 2.0) )
        {
        throw RbException("The dimension between the standard deviations and the correlation coefficients does not match.");
        }
    for (int i=0; i<ce->getValue().size(); i++)
        {
        if ( ce->getValue()[i] < -1.0 || ce->getValue()[i] > 1.0)
            {
            throw RbException("A correlation coefficient is out-of-bounds (-1,+1).");
            }
        }
    RevBayesCore::VarianceCovarianceFunction* f = new RevBayesCore::VarianceCovarianceFunction( sd, ce );
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_varianceCovarianceMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "standardDeviations"  , ModelVector<RealPos>::getClassTypeSpec(), "The vector of standard deviations.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "correlationCoefficients", ModelVector<Real>::getClassTypeSpec(), "The correlation coefficients.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_varianceCovarianceMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_varianceCovarianceMatrix";
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_varianceCovarianceMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_varianceCovarianceMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnVarCovar";
    
    return f_name;
}


const TypeSpec& Func_varianceCovarianceMatrix::getTypeSpec(void) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
