#include "Dist_PhyloMultivariateBrownianREML.h"
#include "RlMatrixRealSymmetric.h"
#include "OptionRule.h"
#include "PhyloMultivariateBrownianProcessREML.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "RlTree.h"

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

/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_PhyloMultivariateBrownianREML::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael R. May" );
    authors.push_back( "Nicolai Vetr" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_PhyloMultivariateBrownianREML::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Multivariate Brownian motion over a phylogeny with variance-covariance matrix rateMatrix.";
    
    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_PhyloMultivariateBrownianREML::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details;
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_PhyloMultivariateBrownianREML::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "\n# generate a tree and variance-covariance matrix.\n";
    example += "psi ~ dnUniformTimeTree(1, [taxon(\"A\"),taxon(\"B\"),taxon(\"C\")])\n";
    example += "Sigma <- diagonalMatrix(5)\n\n";
    example += "# generate the multivariate data.\n";
    example += "x ~ dnPhyloMultivariateBrownianREML(tree=psi, rateMatrix=Sigma)\n\n";
    example += "# print the simulated data.\n";
    example += "x\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_PhyloMultivariateBrownianREML::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    RevBayesCore::RbHelpReference ref1;
    ref1.setCitation("Huelsenbeck JP, B Rannala (2003). Detecting correlation between characters in a comparative analysis with uncertain phylogeny. Evolution, 57(6):1237-1247.");
    references.push_back( ref1 );

    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_PhyloMultivariateBrownianREML::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnPhyloBrownianREML" );
    see_also.push_back( "dnPhyloBrownianMVN" );
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_PhyloMultivariateBrownianREML::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Phylogenetic Multivariate Brownian Motion";
    
    return title;
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

