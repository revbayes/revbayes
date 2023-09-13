#include <math.h>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_multivariateNorm.h"
#include "MemberProcedure.h"
#include "MultivariateNormalDistribution.h"
#include "ModelVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbException.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "RlUtils.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class MatrixReal; }

using namespace RevLanguage;


/**
 * Default constructor.
 * 
 * The default constructor does nothing except allocating the object.
 */
Dist_multivariateNorm::Dist_multivariateNorm(void) : TypedDistribution<ModelVector<Real> >()
{

    
    // member functions
    ArgumentRules* clampAtArgRules = new ArgumentRules();
    clampAtArgRules->push_back( new ArgumentRule( "index", Natural::getClassTypeSpec(), "The index of the value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    clampAtArgRules->push_back( new ArgumentRule( "value", Real::getClassTypeSpec(), "The observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
//    methods.addFunction("clampAt", new DistributionMemberFunction<TimeTree,RealPos>(this, clampAtArgRules   ) );
    methods.addFunction( new MemberProcedure( "clampAt", RlUtils::Void, clampAtArgRules ) );

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
RevBayesCore::MultivariateNormalDistribution* Dist_multivariateNorm::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* m = static_cast<const ModelVector<Real> &>( mean->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* cov = NULL;
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* pre = NULL;
    RevBayesCore::TypedDagNode<double>* sc = static_cast<const RealPos &>( scale->getRevObject() ).getDagNode();
    
    if ( covariance->getRevObject() != RevNullObject::getInstance() && precision->getRevObject() != RevNullObject::getInstance() )
    {
        throw RbException("You can only provide a covariance matrix OR a precision matrix");
    }
    else if ( covariance->getRevObject() != RevNullObject::getInstance() )
    {
        cov = static_cast<const MatrixRealSymmetric &>( covariance->getRevObject() ).getDagNode();
    }
    else if ( precision->getRevObject() != RevNullObject::getInstance() )
    {
        pre = static_cast<const MatrixRealSymmetric &>( precision->getRevObject() ).getDagNode();
    }
    else
    {
        throw RbException("You need to provide a covariance matrix OR a precision matrix");
    }
    
    RevBayesCore::MultivariateNormalDistribution* d = new RevBayesCore::MultivariateNormalDistribution(m, cov, pre, sc);
    
    return d;
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process. 
 */
Dist_multivariateNorm* Dist_multivariateNorm::clone( void ) const
{
    
    return new Dist_multivariateNorm(*this);
}

/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_multivariateNorm::getClassType(void)
{ 
    
    static std::string rev_type = "Dist_multivariateNorm";
    
	return rev_type; 
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_multivariateNorm::getClassTypeSpec(void)
{ 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<ModelVector<Real> >::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_multivariateNorm::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MultivariateNormal";
    
    return d_name;
}


/** 
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the branch rate jump process are:
 * (1) the mean of the distribution.
 * (2) the standard deviation.
 *
 * \return The member rules.
 */
const MemberRules& Dist_multivariateNorm::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set ) 
    {
        dist_member_rules.push_back( new ArgumentRule( "mean"      , ModelVector<Real>::getClassTypeSpec()  , "The vector of mean values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY) );
        dist_member_rules.push_back( new ArgumentRule( "covariance", MatrixRealSymmetric::getClassTypeSpec(), "The variance-covariance matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "precision" , MatrixRealSymmetric::getClassTypeSpec(), "The precision matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "scale"     , RealPos::getClassTypeSpec()            , "The scaling factor of the variance matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_multivariateNorm::getTypeSpec( void ) const
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
void Dist_multivariateNorm::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "mean" ) 
    {
        mean = var;
    }
    else if ( name == "covariance" )
    {
        covariance = var;
    }
    else if ( name == "precision" )
    {
        precision = var;
    }
    else if ( name == "scale" )
    {
        scale = var;
    }
    else 
    {
        Distribution::setConstParameter(name, var);
    }
}
