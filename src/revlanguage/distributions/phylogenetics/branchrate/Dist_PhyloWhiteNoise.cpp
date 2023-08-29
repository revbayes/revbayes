#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_PhyloWhiteNoise.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "PhyloWhiteNoiseProcess.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
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
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

Dist_PhyloWhiteNoise::Dist_PhyloWhiteNoise() : TypedDistribution< ModelVector<RealPos> >()
{
    
}


Dist_PhyloWhiteNoise::~Dist_PhyloWhiteNoise()
{
    
}



Dist_PhyloWhiteNoise* Dist_PhyloWhiteNoise::clone( void ) const
{
    return new Dist_PhyloWhiteNoise(*this);
}


RevBayesCore::PhyloWhiteNoiseProcess* Dist_PhyloWhiteNoise::createDistribution( void ) const
{
    // get the parameters

    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    
    RevBayesCore::TypedDagNode<double>* s  = static_cast<const RealPos &>( sigma->getRevObject() ).getDagNode();
    
    RevBayesCore::PhyloWhiteNoiseProcess* d    = new RevBayesCore::PhyloWhiteNoiseProcess( tau, s );
    
    return d;

}



/* Get Rev type of object */
const std::string& Dist_PhyloWhiteNoise::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloWhiteNoise";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloWhiteNoise::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_PhyloWhiteNoise::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloWhiteNoise";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloWhiteNoise::getParameterRules(void) const
{
    
    static MemberRules dist;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist.push_back( new ArgumentRule( "tree" , TimeTree::getClassTypeSpec(), "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist.push_back( new ArgumentRule( "sigma", RealPos::getClassTypeSpec() , "The standard deviation.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return dist;
}


const TypeSpec& Dist_PhyloWhiteNoise::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}




/** Print value for user */

 void Dist_PhyloWhiteNoise::printValue(std::ostream& o) const
{
    
    o << " whitenoise(";
    
    o << "tau=";
    if ( tree != NULL ) {
        o << tree->getName();
    } else {
        o << "?";
    }
    
    o << "sigma=";
    if ( sigma != NULL ) {
        o << sigma->getName();
    } else {
        o << "?";
    }
    o << ")";
}



/** Set a member variable */
void Dist_PhyloWhiteNoise::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "sigma" )
    {
        sigma = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}

