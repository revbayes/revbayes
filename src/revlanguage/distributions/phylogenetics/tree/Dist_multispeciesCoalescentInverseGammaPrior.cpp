#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Clade.h"
#include "Dist_multispeciesCoalescentInverseGammaPrior.h"
#include "ModelVector.h"
#include "MultispeciesCoalescent.h"
#include "Natural.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "RlClade.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_multispeciesCoalescentInverseGammaPrior::Dist_multispeciesCoalescentInverseGammaPrior() : TypedDistribution<TimeTree>()
{

}

/**
 * Clone the object
 *
 * \return a clone of the object.
 */

Dist_multispeciesCoalescentInverseGammaPrior* Dist_multispeciesCoalescentInverseGammaPrior::clone(void) const
{

    return new Dist_multispeciesCoalescentInverseGammaPrior(*this);

}




/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::MultispeciesCoalescentInverseGammaPrior* Dist_multispeciesCoalescentInverseGammaPrior::createDistribution( void ) const
{

    // Get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* st = static_cast<const TimeTree &>( species_tree->getRevObject() ).getDagNode();
    const std::vector<RevBayesCore::Taxon>      &t  = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    RevBayesCore::MultispeciesCoalescentInverseGammaPrior*   d = new RevBayesCore::MultispeciesCoalescentInverseGammaPrior( st, t );

    RevBayesCore::TypedDagNode<double>* s = static_cast<const RealPos &>( shape->getRevObject() ).getDagNode();
    d->setShape( s );
    RevBayesCore::TypedDagNode<double>* r = static_cast<const RealPos &>( rate->getRevObject() ).getDagNode();
    d->setRate( r );


    d->redrawValue();

    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_multispeciesCoalescentInverseGammaPrior::getClassType(void)
{

    static std::string rev_type = "Dist_multispeciesCoalescentInverseGammaPrior";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_multispeciesCoalescentInverseGammaPrior::getClassTypeSpec(void)
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
std::string Dist_multispeciesCoalescentInverseGammaPrior::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MultiSpeciesCoalescentInverseGamma";

    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Multispecies Coalescent process are:
 * (1) Species tree.
 * (2) Max Population size.
 * (3) Gene name to species name map.
 *
 * \return The member rules.
 */
const MemberRules& Dist_multispeciesCoalescentInverseGammaPrior::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "speciesTree"  , TimeTree::getClassTypeSpec(), "The species tree in which the gene trees evolve.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "shape"        , RealPos::getClassTypeSpec(), "The shape of the inverse gamma prior distribution on the effective population sizes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "rate"         , RealPos::getClassTypeSpec(), "The rate of the inverse gamma prior distribution on the effective population sizes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "taxa"         , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Dist_multispeciesCoalescentInverseGammaPrior::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "speciesTree" )
    {
        species_tree = var;
    }
    else if ( name == "shape" )
    {
        shape = var;
    }
    else if ( name == "rate" )
    {
        rate = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else
    {
        TypedDistribution<TimeTree>::setConstParameter(name, var);
    }

}


/** Get type spec */
const TypeSpec& Dist_multispeciesCoalescentInverseGammaPrior::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_multispeciesCoalescentInverseGammaPrior::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Sebastian Hoehna, Bastien Boussau" );

    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_multispeciesCoalescentInverseGammaPrior::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given effective population sizes. Requires an ultrametric species tree, parameters of an inverse gamma prior on effective population sizes, and taxa with species and individual names.";

    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_multispeciesCoalescentInverseGammaPrior::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "The species tree must be ultrametric.\n";

    details += "Parameters of an inverse gamma prior on effective population sizes must be provided.\n";

    details += "This distribution uses a conjugate prior on effective population sizes.";

    details += " As a consequence, effective population sizes are integrated out and treated as nuisance parameters.\n\n";

    details += "If you are interested in reconstructing ancestral effective population sizes, use dnMultiSpeciesCoalescent.";

    return details;
}


/**
 * Get an executable and instructive example.
 * These examples should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Dist_multispeciesCoalescentInverseGammaPrior::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "# We are going to save the trees we simulate in the folder simulatedTrees:\n";
    example += "dataFolder = \"simulatedTrees/\" \n";
    example += "# Let’s simulate a species tree with 10 taxa, 2 gene trees, 3 alleles per species:\n";
    example += "n_species <- 10\n";
    example += "n_genes <- 2\n";
    example += "n_alleles <- 3\n";
    example += "# we simulate an ultrametric species tree:\n";
    example += "# Species names:\n";
    example += "for (i in 1:n_species) {\n";
    example += "        species[i] <- taxon(taxonName=\"Species_\"+i, speciesName=\"Species_\"+i)\n";
    example += "}\n";
    example += "spTree ~ dnBirthDeath(lambda=0.3, mu=0.2, rootAge=10, rho=1, samplingStrategy=\"uniform\", condition=\"nTaxa\", taxa=species)\n";
    example += "print(spTree)\n";
    example += "# let's pick constant parameters for the inverse gamma distribution:\n";
    example += "alpha <- 3\n";
    example += "beta <- 0.003\n";
    example += "# let's simulate gene trees now:\n";
    example += "# taxa names:\n";
    example += "for (g in 1:n_genes) {\n";
    example += "  for (i in 1:n_species) {\n";
    example += "    for (j in 1:n_alleles) {\n";
    example += "        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName=\"Species_\"+i+\"_\"+j, speciesName=\"Species_\"+i)\n";
    example += "    }\n";
    example += "  }\n";
    example += "  geneTrees[g] ~ dnMultiSpeciesCoalescentInverseGamma(speciesTree=spTree, shape=alpha, rate=beta, taxa=taxons[g])\n";
    example += "  print(geneTrees[g])\n";
    example += "}\n";
    example += "# We can save the species tree and the gene trees: \n";
    example += "write(spTree, filename=dataFolder+\"speciesTree\")\n";
    example += "# Saving the gene trees\n";
    example += "for (i in 1:(n_genes)) {\n";
    example += "  write(geneTrees[i], filename=dataFolder+\"geneTree_\"+i+\".tree\")\n";
    example += "}\n";
    example += "\n";

    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_multispeciesCoalescentInverseGammaPrior::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    RevBayesCore::RbHelpReference ref = RevBayesCore::RbHelpReference();
    ref.setCitation(" Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Jones G.  Journal of Mathematical Biology. 2016.");
    ref.setDoi("DOI: 10.1007/s00285-016-1034-0");
    ref.setUrl("http://www.indriid.com/2016/2016-06-01-STACEY.pdf");

    references.push_back(ref);

    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_multispeciesCoalescentInverseGammaPrior::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnMultiSpeciesCoalescent" );
    see_also.push_back( "dnMultiSpeciesCoalescentUniformPrior" );


    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_multispeciesCoalescentInverseGammaPrior::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Multispecies coalescent Distribution with inverse gamma prior on effective population sizes";

    return title;
}
