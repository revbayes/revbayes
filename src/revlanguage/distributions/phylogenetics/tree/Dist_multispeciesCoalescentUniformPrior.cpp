#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Clade.h"
#include "Dist_multispeciesCoalescentUniformPrior.h"
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
Dist_multispeciesCoalescentUniformPrior::Dist_multispeciesCoalescentUniformPrior() : TypedDistribution<TimeTree>()
{

}

/**
 * Clone the object
 *
 * \return a clone of the object.
 */

Dist_multispeciesCoalescentUniformPrior* Dist_multispeciesCoalescentUniformPrior::clone(void) const
{

    return new Dist_multispeciesCoalescentUniformPrior(*this);

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
RevBayesCore::MultispeciesCoalescentUniformPrior* Dist_multispeciesCoalescentUniformPrior::createDistribution( void ) const
{

    // Get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* st = static_cast<const TimeTree &>( species_tree->getRevObject() ).getDagNode();
    const std::vector<RevBayesCore::Taxon>      &t  = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // get the number of nodes for the tree
    // size_t n_nodes = st->getValue().getNumberOfNodes();


    RevBayesCore::MultispeciesCoalescentUniformPrior*   d = new RevBayesCore::MultispeciesCoalescentUniformPrior( st, t );

    RevBayesCore::TypedDagNode<double>* m = static_cast<const RealPos &>( theta_max->getRevObject() ).getDagNode();
    d->setMaxTheta( m );
    d->redrawValue();

    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_multispeciesCoalescentUniformPrior::getClassType(void)
{

    static std::string rev_type = "Dist_multispeciesCoalescentUniformPrior";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_multispeciesCoalescentUniformPrior::getClassTypeSpec(void)
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
std::string Dist_multispeciesCoalescentUniformPrior::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "MultiSpeciesCoalescentUniformPrior";

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
const MemberRules& Dist_multispeciesCoalescentUniformPrior::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "speciesTree", TimeTree::getClassTypeSpec(), "The species tree in which the gene trees evolve.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "max"    , RealPos::getClassTypeSpec(), "The maximum effective population size.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

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
void Dist_multispeciesCoalescentUniformPrior::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "speciesTree" )
    {
        species_tree = var;
    }
    else if ( name == "max" )
    {
        theta_max = var;
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
const TypeSpec& Dist_multispeciesCoalescentUniformPrior::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



  /**
   * Get the author(s) of this function so they can receive credit (and blame) for it.
   */
  std::vector<std::string> Dist_multispeciesCoalescentUniformPrior::getHelpAuthor(void) const
  {
      // create a vector of authors for this function
      std::vector<std::string> authors;
      authors.push_back( "Sebastian Hoehna, Bastien Boussau" );

      return authors;
  }


  /**
   * Get the (brief) description for this function
   */
  std::string Dist_multispeciesCoalescentUniformPrior::getHelpDescription(void) const
  {
      // create a variable for the description of the function
      std::string description = "Multispecies coalescent distribution describing how gene trees can be generated from within a species tree given effective population sizes. Requires an ultrametric species tree, effective population size(s) (a single real positive or a vector of real positives), and taxa with species and individual names.";

      return description;
  }


  /**
   * Get the more detailed description of the function
   */
  std::string Dist_multispeciesCoalescentUniformPrior::getHelpDetails(void) const
  {
      // create a variable for the description of the function
      std::string details = "";
      details += "The species tree must be ultrametric.\n";

      details += "Effective population sizes can be constant across the species tree, if a single real positive is provided, or branchwise, if a vector is provided. ";

      return details;
  }


  /**
   * Get an executable and instructive example.
   * These examples should help the users to show how this function works but
   * are also used to test if this function still works.
   */
  std::string Dist_multispeciesCoalescentUniformPrior::getHelpExample(void) const
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
      example += "# let's pick a constant effective population size of 50:\n";
      example += "popSize <- 50\n";
      example += "# let's simulate gene trees now:\n";
      example += "# taxa names:\n";
      example += "for (g in 1:n_genes) {\n";
      example += "  for (i in 1:n_species) {\n";
      example += "    for (j in 1:n_alleles) {\n";
      example += "        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName=\"Species_\"+i+\"_\"+j, speciesName=\"Species_\"+i)\n";
      example += "    }\n";
      example += "  }\n";
      example += "  geneTrees[g] ~ dnMultiSpeciesCoalescentUniformPrior(speciesTree=spTree, Ne=popSize, taxa=taxons[g])\n";
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
  std::vector<RevBayesCore::RbHelpReference> Dist_multispeciesCoalescentUniformPrior::getHelpReferences(void) const
  {
      // create an entry for each reference
      std::vector<RevBayesCore::RbHelpReference> references;

      RevBayesCore::RbHelpReference ref = RevBayesCore::RbHelpReference();
      ref.setCitation("Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences From Multiple Loci. Bruce Rannala and Ziheng Yang. GENETICS August 1, 2003 vol. 164 no. 4 1645-1656.");
      ref.setDoi("");
      ref.setUrl("http://www.genetics.org/content/164/4/1645.short");

      references.push_back(ref);

      RevBayesCore::RbHelpReference ref2 = RevBayesCore::RbHelpReference();
      ref2.setCitation("Bayesian Inference of Species Trees from Multilocus Data. Heled and Drummond Mol. Biol Evol. 27 (3): 570-580, 2010.");
      ref2.setDoi("DOI: https://doi.org/10.1093/molbev/msp274");
      ref2.setUrl("https://academic.oup.com/mbe/article/27/3/570/999753/Bayesian-Inference-of-Species-Trees-from");

      references.push_back(ref2);

      return references;
  }


  /**
   * Get the names of similar and suggested other functions
   */
  std::vector<std::string> Dist_multispeciesCoalescentUniformPrior::getHelpSeeAlso(void) const
  {
      // create an entry for each suggested function
      std::vector<std::string> see_also;
      see_also.push_back( "dnMultiSpeciesCoalescent" );
      see_also.push_back( "dnMultiSpeciesCoalescentInverseGamma" );


      return see_also;
  }


  /**
   * Get the title of this help entry
   */
  std::string Dist_multispeciesCoalescentUniformPrior::getHelpTitle(void) const
  {
      // create a title variable
      std::string title = "Multispecies coalescent Distribution";

      return title;
  }
