#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_uniformTimeTree.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "UniformTimeTreeDistribution.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistribution.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "Taxon.h"
#include "Tree.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

using namespace RevLanguage;


Dist_uniformTimeTree::Dist_uniformTimeTree() : TypedDistribution<TimeTree>()
{
    
}


Dist_uniformTimeTree::~Dist_uniformTimeTree()
{

}



Dist_uniformTimeTree* Dist_uniformTimeTree::clone( void ) const
{

    return new Dist_uniformTimeTree( *this );
}


RevBayesCore::UniformTimeTreeDistribution* Dist_uniformTimeTree::createDistribution( void ) const
{

    // Get the parameters
    RevBayesCore::TypedDagNode<double>* a   = static_cast<const RealPos &>( root_age->getRevObject() ).getDagNode();
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    RevBayesCore::UniformTimeTreeDistribution*   d = new RevBayesCore::UniformTimeTreeDistribution( a, t );

    return d;
}



/* Get Rev type of object */
const std::string& Dist_uniformTimeTree::getClassType(void)
{
    
    static std::string rev_type = "Dist_uniformTimeTree";
    
	return rev_type;
}


/* Get class type spec describing type of object. TODO: Check if the correct parent is TypedDistribution or Distribution */
const TypeSpec& Dist_uniformTimeTree::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<TimeTree>::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_uniformTimeTree::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "UniformTimeTree";
    
    return d_name;
}


/* Return member rules */
const MemberRules& Dist_uniformTimeTree::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {

        dist_member_rules.push_back( new ArgumentRule( "rootAge"  , RealPos::getClassTypeSpec()            , "The age of the root.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"     , ModelVector<Taxon>::getClassTypeSpec() , "The taxa used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_uniformTimeTree::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Set a member variable */
void Dist_uniformTimeTree::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "rootAge" )
    {
        root_age = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}

