#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ApproximateTreeLikelihood.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_ApproximateTreeLikelihood.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "Distribution.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RbVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlDistributionMemberFunction.h"
#include "RlMatrixRealSymmetric.h"
#include "RlStochasticNode.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TraceNumeric.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_ApproximateTreeLikelihood::Dist_ApproximateTreeLikelihood() : TypedDistribution<BranchLengthTree>()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_ApproximateTreeLikelihood* Dist_ApproximateTreeLikelihood::clone( void ) const
{

    return new Dist_ApproximateTreeLikelihood(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::ApproximateTreeLikelihood* Dist_ApproximateTreeLikelihood::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::Tree >* tt = static_cast<const TimeTree &>( time_tree->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* br = static_cast<const ModelVector<RealPos> &>( branch_rates->getRevObject() ).getDagNode();
    const RevBayesCore::RbVector<double>* gr = &static_cast<const ModelVector<RealPos> &>( gradients->getRevObject() ).getValue();
    const RevBayesCore::MatrixReal* hess = &static_cast<const MatrixRealSymmetric &>( hessian->getRevObject() ).getDagNode()->getValue();

    const std::string& tr_str = static_cast<const RlString &>( transform->getRevObject() ).getValue();
    RevBayesCore::TRANSFORMATION tr = RevBayesCore::NONE;
    if ( tr_str == "log" ) tr = RevBayesCore::LOG;
    else if ( tr_str == "sqrt" ) tr = RevBayesCore::SQRT;
    else if ( tr_str == "arcsin" ) tr = RevBayesCore::ARCSIN;

    // create the internal distribution object
    RevBayesCore::ApproximateTreeLikelihood* dist = new RevBayesCore::ApproximateTreeLikelihood(tt, br, gr, hess, tr);

    return dist;
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_ApproximateTreeLikelihood::getClassType( void )
{

    static std::string rev_type = "Dist_ApproximateTreeLikelihood";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_ApproximateTreeLikelihood::getClassTypeSpec( void )
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<BranchLengthTree>::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor.
 *
 * \return Rev name of constructor function.
 */
std::string Dist_ApproximateTreeLikelihood::getDistributionFunctionName( void ) const
{

    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "ApproximateTreeLikelihood";
    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the constant-rate birth-death process are:
 * (1) the speciation rate lambda which must be a positive real.
 * (2) the extinction rate mu that must be a positive real.
 * (3) all member rules specified by BirthDeathProcess.
 *
 * \return The member rules.
 */
const MemberRules& Dist_ApproximateTreeLikelihood::getParameterRules(void) const
{

    static MemberRules member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        member_rules.push_back( new ArgumentRule( "timeTree", TimeTree::getClassTypeSpec(), "The time tree", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "branchRates", ModelVector<RealPos>::getClassTypeSpec(), "The vector of branch rates.",   ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "gradientVector", ModelVector<RealPos>::getClassTypeSpec(), "The gradient vector (vector of first derivatives of the log-likelihood function).", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        member_rules.push_back( new ArgumentRule( "hessianMatrix", MatrixRealSymmetric::getClassTypeSpec(), "The Hessian matrix (matrix of second derivatives of the log-likelihood functon).", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<std::string> options;
        options.push_back( "none" );
        options.push_back( "log" );
        options.push_back( "sqrt" );
        options.push_back( "arcsin" );
        
        member_rules.push_back( new OptionRule( "branchLengthTransform", new RlString("none"), options, "What transformation should we apply to branch lengths in expected substitutions per character?" ) );

        rules_set = true;
    }
    
    return member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_ApproximateTreeLikelihood::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_ApproximateTreeLikelihood::setConstParameter(const std::string& name, const RevPtr<const RevVariable>& var)
{

    if ( name == "timeTree" )
    {
        time_tree = var;
    }
    else if ( name == "branchRates" )
    {
        branch_rates = var;
    }
    else if ( name == "gradientVector" )
    {
        gradients = var;
    }
    else if ( name == "hessianMatrix" )
    {
        hessian = var;
    }
    else if ( name == "branchLengthTransform" )
    {
        transform = var;
    }
    else
    {
        TypedDistribution<BranchLengthTree>::setConstParameter(name, var);
    }
}
