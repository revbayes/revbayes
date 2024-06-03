#include <math.h>
#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_SFSDiffusionApproximation.h"
#include "SFSDiffusionApproximationDistribution.h"
#include "RbVector.h"
#include "RlSimplex.h"
#include "ConstantNode.h"
#include "DagMemberFunction.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DistributionMemberFunction.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistributionMemberFunction.h"
#include "RlStochasticNode.h"
#include "RlString.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;

Dist_SFSDiffusionApproximation::Dist_SFSDiffusionApproximation(void) : TypedDistribution< ModelVector<RealPos> >()
{
    
}


Dist_SFSDiffusionApproximation::~Dist_SFSDiffusionApproximation(void)
{
    
}



Dist_SFSDiffusionApproximation* Dist_SFSDiffusionApproximation::clone( void ) const
{
    
    return new Dist_SFSDiffusionApproximation(*this);
}


RevBayesCore::SFSDiffusionApproximationDistribution* Dist_SFSDiffusionApproximation::createDistribution( void ) const
{
    
    // get the parameters
    long                                                            n_sites = static_cast<const Natural              &>( num_sites->getRevObject() ).getValue();
    long                                                            n_ind   = static_cast<const Natural              &>( num_individuals->getRevObject() ).getValue();
    long                                                            f       = static_cast<const RlBoolean            &>( folded->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*   th      = static_cast<const ModelVector<RealPos> &>( theta->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*   ls      = static_cast<const ModelVector<RealPos> &>( lengths->getRevObject() ).getDagNode();

    const std::string&                                              c       = static_cast<const RlString             &>( coding->getRevObject() ).getValue();
    RevBayesCore::SFSDiffusionApproximationDistribution::CODING cd = RevBayesCore::SFSDiffusionApproximationDistribution::ALL;
    if ( c == "no-monomorphic" )
    {
        cd = RevBayesCore::SFSDiffusionApproximationDistribution::NO_MONOMORPHIC;
    }
    else if ( c == "no-singletons" )
    {
        cd = RevBayesCore::SFSDiffusionApproximationDistribution::NO_SINGLETONS;
    }

    RevBayesCore::SFSDiffusionApproximationDistribution*            d       = new RevBayesCore::SFSDiffusionApproximationDistribution( th, ls, n_sites, n_ind, f, cd );
    
    return d;
}



/* Get Rev type of object */
const std::string& Dist_SFSDiffusionApproximation::getClassType(void)
{
    static std::string rev_type = "Dist_SFSDiffusionApproximation";
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_SFSDiffusionApproximation::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<RealPos> >::getClassTypeSpec() ) );
    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_SFSDiffusionApproximation::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "SFSDiffusionApproximation";
    
    return d_name;
}


MethodTable Dist_SFSDiffusionApproximation::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution< ModelVector<RealPos> >::getDistributionMethods();
    
    ArgumentRules* esfs_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_SFSDiffusionApproximation, ModelVector<RealPos> >( "getExpectedAlleleFrequencies", variable, esfs_arg_rules, true ) );

    return methods;
}



/** Return member rules (no members) */
const MemberRules& Dist_SFSDiffusionApproximation::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "theta", ModelObject<RealPos>::getClassTypeSpec(), "The theta values with theta=4*Ne*mu.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "lengths", ModelObject<RealPos>::getClassTypeSpec(), "The epoch lengths for all but the most ancient epoch. Time is in units of 2*Na generations with Na being the most ancient population size.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "numSites", Natural::getClassTypeSpec(), "The number of sites in the SFS.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "numIndividuals", Natural::getClassTypeSpec(), "The number of individuals in the (unfolded) SFS.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "folded", RlBoolean::getClassTypeSpec(), "Is the site frequency spectrum folded?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(false) ) );

        std::vector<std::string> coding_options;
        coding_options.push_back( "all" );
        coding_options.push_back( "no-monomorphic" );
        coding_options.push_back( "no-singletons" );
        dist_member_rules.push_back( new OptionRule( "coding", new RlString("all"), coding_options, "The assumption which allele frequencies are included." ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_SFSDiffusionApproximation::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_SFSDiffusionApproximation::printValue(std::ostream& o) const
{
    o << "SFSDiffusionApproximation(theta=";
    if ( theta != NULL )
        o << theta->getName();
    else
        o << "?";
    o << ")";
}


/** Set a member variable */
void Dist_SFSDiffusionApproximation::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "theta" )
    {
        theta = var;
    }
    else if ( name == "lengths" )
    {
        lengths = var;
    }
    else if ( name == "numSites" )
    {
        num_sites = var;
    }
    else if ( name == "numIndividuals" )
    {
        num_individuals = var;
    }
    else if ( name == "folded" )
    {
        folded = var;
    }
    else if ( name == "coding" )
    {
        coding = var;
    }
    else
    {
        TypedDistribution< ModelVector<RealPos> >::setConstParameter(name, var);
    }
}
