#include <math.h>
#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_EpochCoalRateMatrixDemography.h"
#include "EpochCoalRateMatrixDemography.h"
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
#include "RlSimplex.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class Simplex; }

using namespace RevLanguage;

Dist_EpochCoalRateMatrixDemography::Dist_EpochCoalRateMatrixDemography(void) : TypedDistribution< ModelVector<RealPos> >()
{
    
}


Dist_EpochCoalRateMatrixDemography::~Dist_EpochCoalRateMatrixDemography(void)
{
    
}



Dist_EpochCoalRateMatrixDemography* Dist_EpochCoalRateMatrixDemography::clone( void ) const
{
    
    return new Dist_EpochCoalRateMatrixDemography(*this);
}


RevBayesCore::EpochCoalRateMatrixDemography* Dist_EpochCoalRateMatrixDemography::createDistribution( void ) const
{
    
    // get the parameters
    long                                                            n_sites = static_cast<const Natural              &>( num_sites->getRevObject() ).getValue();
    long                                                            n_ind   = static_cast<const Natural              &>( num_individuals->getRevObject() ).getValue();
    long                                                            f       = static_cast<const RlBoolean            &>( folded->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*   ne      = static_cast<const ModelVector<RealPos> &>( Ne->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*   et      = static_cast<const ModelVector<RealPos> &>( times->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >*   m       = static_cast<const ModelVector<RealPos> &>( mu->getRevObject() ).getDagNode();
    RevBayesCore::TypedDagNode< RevBayesCore::Simplex >*            asfs    = static_cast<const Simplex              &>( ancestral_sfs->getRevObject() ).getDagNode();

    const std::string&                                              c       = static_cast<const RlString             &>( coding->getRevObject() ).getValue();
    RevBayesCore::EpochCoalRateMatrixDemography::CODING cd = RevBayesCore::EpochCoalRateMatrixDemography::ALL;
    if ( c == "no-monomorphic" )
    {
        cd = RevBayesCore::EpochCoalRateMatrixDemography::NO_MONOMORPHIC;
    }
    else if ( c == "no-singletons" )
    {
        cd = RevBayesCore::EpochCoalRateMatrixDemography::NO_SINGLETONS;
    }

    RevBayesCore::EpochCoalRateMatrixDemography* d = new RevBayesCore::EpochCoalRateMatrixDemography( ne, et, m, asfs, n_sites, n_ind, f, cd );
       
    return d;
}



/* Get Rev type of object */
const std::string& Dist_EpochCoalRateMatrixDemography::getClassType(void)
{
    
    static std::string rev_type = "Dist_EpochCoalRateMatrixDemography";
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_EpochCoalRateMatrixDemography::getClassTypeSpec(void)
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
std::string Dist_EpochCoalRateMatrixDemography::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "EpochCoalRateMatrixDemography";
    
    return d_name;
}


MethodTable Dist_EpochCoalRateMatrixDemography::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution< ModelVector<RealPos> >::getDistributionMethods();

    ArgumentRules* esfs_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_EpochCoalRateMatrixDemography, ModelVector<RealPos> >( "getExpectedAlleleFrequencies", variable, esfs_arg_rules, true ) );

    return methods;
}



/** Return member rules (no members) */
const MemberRules& Dist_EpochCoalRateMatrixDemography::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "Ne", ModelObject<RealPos>::getClassTypeSpec(), "The effective population size per epoch.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "times", ModelObject<RealPos>::getClassTypeSpec(), "The times when the epochs end.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "ancestralSFS", Simplex::getClassTypeSpec(), "The ancestral probabilities of site frequencies.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "mu", ModelObject<RealPos>::getClassTypeSpec(), "The mutation rates.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "numSites", Natural::getClassTypeSpec(), "The number of sites in the SFS.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "numIndividuals", Natural::getClassTypeSpec(), "The number of individuals in (unfolded) the SFS.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        dist_member_rules.push_back( new ArgumentRule( "folded", RlBoolean::getClassTypeSpec(), "Is the site frequency folded.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlBoolean(false) ) );

        std::vector<std::string> coding_options;
        coding_options.push_back( "all" );
        coding_options.push_back( "no-monomorphic" );
        coding_options.push_back( "no-singletons" );
        dist_member_rules.push_back( new OptionRule( "coding", new RlString("all"), coding_options, "The assumption which allele frequencies are included." ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_EpochCoalRateMatrixDemography::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


/** Print value for user */
void Dist_EpochCoalRateMatrixDemography::printValue(std::ostream& o) const
{
    
    o << "EpochCoalRateMatrixDemography(Ne=";
    if ( Ne != NULL )
        o << Ne->getName();
    else
        o << "?";
    o << ")";
}


/** Set a member variable */
void Dist_EpochCoalRateMatrixDemography::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "Ne" )
    {
        Ne = var;
    }
    else if ( name == "times" )
    {
        times = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "ancestralSFS" )
    {
        ancestral_sfs = var;
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
