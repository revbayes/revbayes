#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ScaledDirichletDistribution.h"
#include "Dist_scaledDirichlet.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
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
#include "Simplex.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

Dist_scaledDirichlet::Dist_scaledDirichlet() : TypedDistribution<Simplex>()
{
    
}


Dist_scaledDirichlet::~Dist_scaledDirichlet()
{
    
}



Dist_scaledDirichlet* Dist_scaledDirichlet::clone( void ) const
{
    return new Dist_scaledDirichlet(*this);
}


RevBayesCore::ScaledDirichletDistribution* Dist_scaledDirichlet::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<RealPos> &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* b = static_cast<const ModelVector<RealPos> &>( beta->getRevObject() ).getDagNode();
    RevBayesCore::ScaledDirichletDistribution* d              = new RevBayesCore::ScaledDirichletDistribution( a, b );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_scaledDirichlet::getClassType(void)
{
    
    static std::string rev_type = "Dist_scaledDirichlet";
    
	return rev_type; 
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_scaledDirichlet::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Simplex>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_scaledDirichlet::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "scaledDirichlet";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_scaledDirichlet::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "alpha", ModelVector<RealPos>::getClassTypeSpec(), "The concentration parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "beta", ModelVector<RealPos>::getClassTypeSpec(), "The rate parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


const TypeSpec& Dist_scaledDirichlet::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_scaledDirichlet::printValue(std::ostream& o) const
{
    
    o << " dirichlet(alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ",beta=";
    if ( beta != NULL )
    {
        o << beta->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
}


/** Set a member variable */
void Dist_scaledDirichlet::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "alpha" ) 
    {
        alpha = var;
    }
    else if ( name == "beta" )
    {
        beta = var;
    }

    else 
    {
        Distribution::setConstParameter(name, var);
    }
}
