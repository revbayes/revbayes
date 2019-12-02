//
//  Func_varianceCovarianceMatrix.cpp
//  RevBayesCore
//
//  Created by John Huelsenbeck on 11/17/14.
//

#include "Func_decomposedVarianceCovarianceMatrix.h"
#include "MatrixReal.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "TypedDagNode.h"
#include "DecomposedVarianceCovarianceFunction.h"

using namespace RevLanguage;


/** default constructor */
Func_decomposedVarianceCovarianceMatrix::Func_decomposedVarianceCovarianceMatrix(void) : TypedFunction<MatrixRealSymmetric>()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_decomposedVarianceCovarianceMatrix* Func_decomposedVarianceCovarianceMatrix::clone(void) const {
    
    return new Func_decomposedVarianceCovarianceMatrix(*this);
}


RevBayesCore::TypedFunction< RevBayesCore::MatrixReal >* Func_decomposedVarianceCovarianceMatrix::createFunction(void) const {

    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sd = static_cast<const ModelVector<RealPos>& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* r = static_cast<const MatrixReal& >( this->args[1].getVariable()->getRevObject() ).getDagNode();
    if ( sd->getValue().size() != r->getValue().getNumberOfColumns() )
        {
        throw RbException("The dimension between the standard deviations and the correlation coefficients does not match.");
        }
    for (int i=0; i<r->getValue().getNumberOfRows(); i++)
        {
        for (int j=i+1; j<r->getValue().getNumberOfColumns(); j++)
            {
            if ( r->getValue()[i][j] < -1.0 || r->getValue()[i][j] > 1.0)
                {
                throw RbException("A correlation coefficient is out-of-bounds (-1,+1).");
                }
            }
        }
    RevBayesCore::DecomposedVarianceCovarianceFunction* f = new RevBayesCore::DecomposedVarianceCovarianceFunction( sd, r );
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_decomposedVarianceCovarianceMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "standardDeviations"  , ModelVector<RealPos>::getClassTypeSpec(), "A vector of standard deviations.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "correlationCoefficients", MatrixReal::getClassTypeSpec(), "A matrix of correlation coefficients.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    return argumentRules;
}


const std::string& Func_decomposedVarianceCovarianceMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_decomposedVarianceCovarianceMatrix";
	return rev_type;
}


/* Get class type spec describing type of object */
const TypeSpec& Func_decomposedVarianceCovarianceMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_decomposedVarianceCovarianceMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDecompVarCovar";
    
    return f_name;
}


const TypeSpec& Func_decomposedVarianceCovarianceMatrix::getTypeSpec(void) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
