#include <math.h>
#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "SortedDirichletDistribution.h"
#include "Dist_sortedDirichlet.h"
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

Dist_sortedDirichlet::Dist_sortedDirichlet() : TypedDistribution<Simplex>()
{

}


Dist_sortedDirichlet::~Dist_sortedDirichlet()
{

}



Dist_sortedDirichlet* Dist_sortedDirichlet::clone( void ) const
{
    return new Dist_sortedDirichlet(*this);
}


RevBayesCore::SortedDirichletDistribution* Dist_sortedDirichlet::createDistribution( void ) const
{

    // get the parameters
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<RealPos> &>( alpha->getRevObject() ).getDagNode();
    RevBayesCore::SortedDirichletDistribution* d              = new RevBayesCore::SortedDirichletDistribution( a );

    return d;
}



/* Get Rev type of object */
const std::string& Dist_sortedDirichlet::getClassType(void)
{

    static std::string rev_type = "Dist_sortedDirichlet";

	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_sortedDirichlet::getClassTypeSpec(void)
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
std::string Dist_sortedDirichlet::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "sortedDirichlet";

    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_sortedDirichlet::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "alpha", ModelVector<RealPos>::getClassTypeSpec(), "The concentration parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return memberRules;
}


const TypeSpec& Dist_sortedDirichlet::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** Print value for user */
void Dist_sortedDirichlet::printValue(std::ostream& o) const
{

    o << " sortedDirichlet(alpha=";
    if ( alpha != NULL )
    {
        o << alpha->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";

}


/** Set a member variable */
void Dist_sortedDirichlet::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "alpha" )
    {
        alpha = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
}
