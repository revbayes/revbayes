#include "Dist_PhyloBrownianMultiSampleREML.h"

#include <cstddef>
#include <ostream>

#include "PhyloBrownianProcessMultiSampleREML.h"
#include "RlTaxon.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RbException.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlDistribution.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Taxon; }

using namespace RevLanguage;

Dist_PhyloBrownianMultiSampleREML::Dist_PhyloBrownianMultiSampleREML() : TypedDistribution< ContinuousCharacterData >()
{
    
}


Dist_PhyloBrownianMultiSampleREML::~Dist_PhyloBrownianMultiSampleREML()
{
    
}



Dist_PhyloBrownianMultiSampleREML* Dist_PhyloBrownianMultiSampleREML::clone( void ) const
{
    
    return new Dist_PhyloBrownianMultiSampleREML(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloBrownianMultiSampleREML::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    size_t n = size_t( static_cast<const Natural &>( n_sites->getRevObject() ).getValue() );
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* var = static_cast<const ModelVector<RealPos> &>( within_species_variances->getRevObject() ).getDagNode();
    const std::vector<RevBayesCore::Taxon>      &t  = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    
    RevBayesCore::PhyloBrownianProcessMultiSampleREML *dist = new RevBayesCore::PhyloBrownianProcessMultiSampleREML(tau, var, t, n);
    
    // set the clock rates
    if ( branch_rates->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* br = static_cast<const ModelVector<RealPos> &>( branch_rates->getRevObject() ).getDagNode();
        
        // sanity check
        size_t n_rates = br->getValue().size();
        if ( (n_nodes-1) != n_rates )
        {
            throw RbException( "The number of clock rates does not match the number of branches" );
        }
        
        dist->setBranchRate( br );
        
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* br = static_cast<const RealPos &>( branch_rates->getRevObject() ).getDagNode();
        
        dist->setBranchRate( br );
    }
    
    // set the clock rates
    if ( site_rates->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* sr = static_cast<const ModelVector<RealPos> &>( site_rates->getRevObject() ).getDagNode();
        dist->setSiteRate( sr );
    }
    else
    {
        RevBayesCore::TypedDagNode< double >* sr = static_cast<const RealPos &>( site_rates->getRevObject() ).getDagNode();
        dist->setSiteRate( sr );
    }
    
    return dist;
}



/* Get Rev type of object */
const std::string& Dist_PhyloBrownianMultiSampleREML::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloBrownianMultiSampleREML";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloBrownianMultiSampleREML::getClassTypeSpec(void)
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
std::string Dist_PhyloBrownianMultiSampleREML::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloBrownianMultiSampleREML";
    
    return d_name;
}


/** Return member rules (no members) */
const MemberRules& Dist_PhyloBrownianMultiSampleREML::getParameterRules(void) const
{
    
    static MemberRules dist_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        dist_member_rules.push_back( new ArgumentRule( "tree" , Tree::getClassTypeSpec(), "The tree along which the process evolves.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> branchRateTypes;
        branchRateTypes.push_back( RealPos::getClassTypeSpec() );
        branchRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "branchRates" , branchRateTypes, "The per branch rate-multiplier(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );
        
        dist_member_rules.push_back( new ArgumentRule( "withinSpeciesVariances" , ModelVector<RealPos>::getClassTypeSpec(), "The per species within-species variances.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        std::vector<TypeSpec> siteRateTypes;
        siteRateTypes.push_back( RealPos::getClassTypeSpec() );
        siteRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        RealPos *defaultSiteRates = new RealPos(1.0);
        dist_member_rules.push_back( new ArgumentRule( "siteRates" , siteRateTypes, "The per site rate-multiplier(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultSiteRates ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        dist_member_rules.push_back( new ArgumentRule( "nSites"         ,  Natural::getClassTypeSpec(), "The number of sites used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10) ) );
        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_PhyloBrownianMultiSampleREML::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_PhyloBrownianMultiSampleREML::printValue(std::ostream& o) const
{
    
    o << "PhyloBrownianProcess(tree=";
    if ( tree != NULL )
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ", branchRates=";
    if ( branch_rates != NULL )
    {
        o << branch_rates->getName();
    }
    else
    {
        o << "?";
    }
    o << ", siteRates=";
    if ( site_rates != NULL )
    {
        o << site_rates->getName();
    }
    else
    {
        o << "?";
    }
    o << ", nSites=";
    if ( n_sites != NULL )
    {
        o << n_sites->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
    
}


/** Set a member variable */
void Dist_PhyloBrownianMultiSampleREML::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "branchRates" )
    {
        branch_rates = var;
    }
    else if ( name == "siteRates" )
    {
        site_rates = var;
    }
    else if ( name == "withinSpeciesVariances" )
    {
        within_species_variances = var;
    }
    else if ( name == "nSites" )
    {
        n_sites = var;
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

