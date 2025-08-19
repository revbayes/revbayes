#include "Dist_PhyloMultivariateBrownianREML.h"

#include <cstddef>
#include <ostream>

#include "RlMatrixRealSymmetric.h"
#include "PhyloMultivariateBrownianProcessREML.h"
#include "RlTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RbException.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RealPos.h"
#include "RlDistribution.h"
#include "RlMatrixReal.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TypeSpec.h"

using namespace RevLanguage;

Dist_PhyloMultivariateBrownianREML::Dist_PhyloMultivariateBrownianREML() : TypedDistribution< ContinuousCharacterData >()
{
    
}


Dist_PhyloMultivariateBrownianREML::~Dist_PhyloMultivariateBrownianREML()
{
    
}



Dist_PhyloMultivariateBrownianREML* Dist_PhyloMultivariateBrownianREML::clone( void ) const
{
    
    return new Dist_PhyloMultivariateBrownianREML(*this);
}


RevBayesCore::TypedDistribution< RevBayesCore::ContinuousCharacterData >* Dist_PhyloMultivariateBrownianREML::createDistribution( void ) const
{
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
//    size_t n = size_t( static_cast<const Natural &>( nSites->getRevObject() ).getValue() );
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* vcv = static_cast<const MatrixRealSymmetric&>( rate_matrix->getRevObject() ).getDagNode();

    size_t n = vcv->getValue().getDim();
    
    RevBayesCore::PhyloMultivariateBrownianProcessREML *dist = new RevBayesCore::PhyloMultivariateBrownianProcessREML(tau, vcv, n);

    // set the clock rates
    if ( branchRates->getRevObject().isType( ModelVector<RealPos>::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* br = static_cast<const ModelVector<RealPos> &>( branchRates->getRevObject() ).getDagNode();
        
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
        RevBayesCore::TypedDagNode< double >* br = static_cast<const RealPos &>( branchRates->getRevObject() ).getDagNode();
        
        dist->setBranchRate( br );
    }
    
    return dist;
}



/* Get Rev type of object */
const std::string& Dist_PhyloMultivariateBrownianREML::getClassType(void)
{
    
    static std::string rev_type = "Dist_PhyloMultivariateBrownianREML";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_PhyloMultivariateBrownianREML::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Distribution::getClassTypeSpec() ) );
    
    return rev_type_spec;
}

/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_PhyloMultivariateBrownianREML::getDistributionFunctionAliases( void ) const
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
std::string Dist_PhyloMultivariateBrownianREML::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "PhyloMultivariateBrownianREML";
    
    return d_name;
}



/** Return member rules (no members) */
const MemberRules& Dist_PhyloMultivariateBrownianREML::getParameterRules(void) const
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
        
        // std::vector<TypeSpec> siteRateTypes;
        // siteRateTypes.push_back( RealPos::getClassTypeSpec() );
        // siteRateTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
        // RealPos *defaultSiteRates = new RealPos(1.0);
        // dist_member_rules.push_back( new ArgumentRule( "siteRates" , siteRateTypes, "The per site rate-multiplier(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, defaultSiteRates ) );
        
        std::vector<TypeSpec> rateMatrixTypes;
        rateMatrixTypes.push_back( MatrixRealSymmetric::getClassTypeSpec() );
        rateMatrixTypes.push_back( MatrixReal::getClassTypeSpec() );
        dist_member_rules.push_back( new ArgumentRule( "rateMatrix", rateMatrixTypes, "The variance-covariance matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
//        dist_member_rules.push_back( new ArgumentRule( "nSites"         ,  Natural::getClassTypeSpec(), "The number of sites used for simulation.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(10) ) );

        
        rules_set = true;
    }
    
    return dist_member_rules;
}


const TypeSpec& Dist_PhyloMultivariateBrownianREML::getTypeSpec( void ) const
{
    
    static TypeSpec ts = getClassTypeSpec();
    
    return ts;
}


/** Print value for user */
void Dist_PhyloMultivariateBrownianREML::printValue(std::ostream& o) const
{
    
    o << "PhyloMultivariateBrownianProcess(tree=";
    if ( tree != NULL )
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ", branchRates=";
    if ( branchRates != NULL )
    {
        o << branchRates->getName();
    }
    else
    {
        o << "?";
    }
    // o << ", siteRates=";
    // if ( site_rates != NULL )
    // {
    //     o << site_rates->getName();
    // }
    // else
    //{
    //     o << "?";
    // }
//    o << ", nSites=";
//    if ( nSites != NULL )
//    {
//        o << nSites->getName();
//    }
//    else
//    {
//        o << "?";
//    }
    o << ")";
    
}


/** Set a member variable */
void Dist_PhyloMultivariateBrownianREML::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "branchRates" )
    {
        branchRates = var;
    }
    // else if ( name == "siteRates" )
    // {
    //    site_rates = var;
    // }
//    else if ( name == "nSites" )
//    {
//        nSites = var;
//    }
    else if ( name == "rateMatrix" )
    {
        rate_matrix = var;
    }
    else
    {
        Distribution::setConstParameter(name, var);
    }
    
}

