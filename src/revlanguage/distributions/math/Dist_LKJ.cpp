#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Dist_LKJ.h"
#include "Natural.h"
#include "RealPos.h"
#include "RlMatrixRealSymmetric.h"
#include "StochasticNode.h"
#include "LKJDistribution.h"
#include "ArgumentRules.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;

Dist_LKJ::Dist_LKJ() : TypedDistribution<MatrixRealSymmetric>()
{
    
}


Dist_LKJ::~Dist_LKJ()
{
    
}



Dist_LKJ* Dist_LKJ::clone( void ) const
{
    return new Dist_LKJ(*this);
}


RevBayesCore::LKJDistribution* Dist_LKJ::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<double>* et = static_cast<const RealPos&>( eta->getRevObject() ).getDagNode();
    size_t d  = size_t( static_cast<const Natural &>( dim->getRevObject() ).getValue() );
    
    RevBayesCore::LKJDistribution* f = new RevBayesCore::LKJDistribution( et, d );

    return f;

}



/* Get Rev type of object */
const std::string& Dist_LKJ::getClassType(void)
{
    static std::string rev_type = "Dist_LKJ";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_LKJ::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return rev_type_spec;
}

/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_LKJ::getDistributionFunctionAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    
    return a_names;
}



/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_LKJ::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "LKJ";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_LKJ::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dist_member_rules.push_back( new ArgumentRule( "eta", RealPos::getClassTypeSpec(), "The parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "dim", Natural::getClassTypeSpec(), "The dimensions of the correlation matrix.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_LKJ::getTypeSpec( void ) const {
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_LKJ::printValue(std::ostream& o) const {
    
    o << " LKJ(";
    o << ")";
}


/** Set a member variable */
void Dist_LKJ::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "eta" )
    {
        eta = var;
    }
    else if ( name == "dim" )
    {
        dim = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
