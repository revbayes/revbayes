#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_SpeciesSubtreeScale.h"
#include "SpeciesSubtreeScaleProposal.h"
#include "RbException.h"
#include "RealPos.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_SpeciesSubtreeScale::Move_SpeciesSubtreeScale() : Move()
{

    // add method for call "addGeneTreeVariable" as a function
    ArgumentRules* addGeneTreeArgRules = new ArgumentRules();
    addGeneTreeArgRules->push_back( new ArgumentRule( "geneTree" , TimeTree::getClassTypeSpec(), "A gene tree to scale.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    methods.addFunction( new MemberProcedure( "addGeneTreeVariable", RlUtils::Void, addGeneTreeArgRules) );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_SpeciesSubtreeScale* Move_SpeciesSubtreeScale::clone(void) const
{

    return new Move_SpeciesSubtreeScale(*this);
}


/**
 * Create a new internal move object.
 *
 * This function simply dynamically allocates a new internal move object that is
 * associated with the variable (DAG-node). The internal move object is created by calling its
 * constructor and passing the move-parameters (the variable and other parameters) as arguments of the
 * constructor. The move constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
void Move_SpeciesSubtreeScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( speciesTree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *st = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );


    RevBayesCore::Proposal *p = new RevBayesCore::SpeciesSubtreeScaleProposal(st);
    value = new RevBayesCore::MetropolisHastingsMove(p,w);

}


RevPtr<RevVariable> Move_SpeciesSubtreeScale::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{

    if ( name == "addGeneTreeVariable" )
    {
        found = true;

        RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( args[0].getVariable()->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::Tree> *gt = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::SpeciesSubtreeScaleProposal &p = static_cast<RevBayesCore::SpeciesSubtreeScaleProposal&>( m->getProposal() );
        p.addGeneTree( gt );

        return NULL;
    }
    else if ( name == "removeGeneTreeVariable" )
    {
        found = true;

        RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( args[0].getVariable()->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::Tree> *gt = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::SpeciesSubtreeScaleProposal &p = static_cast<RevBayesCore::SpeciesSubtreeScaleProposal&>( m->getProposal() );
        p.removeGeneTree( gt );

        return NULL;
    }

    return Move::executeMethod( name, args, found );
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_SpeciesSubtreeScale::getClassType(void)
{

    static std::string rev_type = "Move_SpeciesSubtreeScale";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_SpeciesSubtreeScale::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SpeciesSubtreeScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SpeciesSubtreeScale";

    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be a time-tree.
 *
 * \return The member rules.
 */
const MemberRules& Move_SpeciesSubtreeScale::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "speciesTree", TimeTree::getClassTypeSpec(), "The species variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC    ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_SpeciesSubtreeScale::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/**
 * Print the value for the user.
 */
void Move_SpeciesSubtreeScale::printValue(std::ostream &o) const
{

    o << "SpeciesSubtreeScale(";
    if (speciesTree != NULL)
    {
        o << speciesTree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";

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
void Move_SpeciesSubtreeScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "speciesTree" )
    {
        speciesTree = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }

}



  /**
   * Get the author(s) of this function so they can receive credit (and blame) for it.
   */
  std::vector<std::string> Move_SpeciesSubtreeScale::getHelpAuthor(void) const
  {
      // create a vector of authors for this function
      std::vector<std::string> authors;
      authors.push_back( "Sebastian Hoehna, Bastien Boussau" );

      return authors;
  }


  /**
   * Get the (brief) description for this function
   */
  std::string Move_SpeciesSubtreeScale::getHelpDescription(void) const
  {
      // create a variable for the description of the function
      std::string description =  "Makes a subtree scale move both in the species tree and in the gene trees that contain nodes of the relevant populations. Tree topologies are not altered." ;

      return description;
  }


  /**
   * Get the more detailed description of the function
   */
  std::string Move_SpeciesSubtreeScale::getHelpDetails(void) const
  {
      // create a variable for the description of the function
      std::string details = "";
      details += "The species tree must be ultrametric.\n\n";

      details += "All the gene trees that evolved along the species tree according to some form of multispecies coalescent must be added to the move using the addGeneTreeVariable method.\n\n";

      details += "This move jointly performs a subtree scale move (a whole subtree is scaled up or down, keeping the topology fixed) on the species tree and on gene trees, all of which must be ultrametric.\n\n";

      details += "How this works: we pick a random node which is not the root.\nThen, we uniformly pick an age between the parent and the oldest sampled descendant.\nThe picked subtree is then scaled to this new age.\nAll gene-trees that are present in the population will be scaled accordingly.";

      return details;
  }


  /**
   * Get an executable and instructive example.
   * These examples should help the users to show how this function works but
   * are also used to test if this function still works.
   */
  std::string Move_SpeciesSubtreeScale::getHelpExample(void) const
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
      example += "  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])\n";
      example += "  print(geneTrees[g])\n";
      example += "}\n";
      example += "# We can save the species tree and the gene trees: \n";
      example += "write(spTree, filename=dataFolder+\"speciesTree\")\n";
      example += "# Saving the gene trees\n";
      example += "for (i in 1:(n_genes)) {\n";
      example += "  write(geneTrees[i], filename=dataFolder+\"geneTree_\"+i+\".tree\")\n";
      example += "}\n";
      example += "# set my move index\n";
      example += "mi = 0\n";
      example += "move_species_subtree_scale = mvSpeciesSubtreeScale( speciesTree=spTree, weight=5 )\n";
      example += "for (i in 1:n_genes) {\n";
      example += "   move_species_subtree_scale.addGeneTreeVariable( geneTrees[i] )\n";
      example += "}\n";
      example += "moves[++mi] = move_species_subtree_scale\n";
      example += "# We get a handle on our model.\n";
      example += "# We can use any node of our model as a handle, here we choose to use the topology.\n";
      example += "mymodel = model(spTree)\n";
      example += "# Monitors to check the progression of the program\n";
      example += "monitors[1] = mnScreen(printgen=10, spTree)\n";
      example += "# Here we use a plain MCMC. You could also set nruns=2 for a replicated analysis\n";
      example += "# or use mcmcmc with heated chains.\n";
      example += "mymcmc = mcmc(mymodel, monitors, moves, nruns=4)\n";
      example += "mymcmc.run(generations=1000)\n";
      example += "mymcmc.operatorSummary()\n";

      return example;
  }


  /**
   * Get some references/citations for this function
   *
   */
  std::vector<RevBayesCore::RbHelpReference> Move_SpeciesSubtreeScale::getHelpReferences(void) const
  {
      // create an entry for each reference
      std::vector<RevBayesCore::RbHelpReference> references;
      RevBayesCore::RbHelpReference ref = RevBayesCore::RbHelpReference();
      ref.setCitation("Guided Tree Topology Proposals for Bayesian Phylogenetic Inference. Sebastian Höhna, Alexei J. Drummond. Syst Biol (2012) 61 (1): 1-11.");
      ref.setDoi("https://doi.org/10.1093/sysbio/syr074");
      ref.setUrl("https://academic.oup.com/sysbio/article-lookup/doi/10.1093/sysbio/syr074");

      references.push_back(ref);

      RevBayesCore::RbHelpReference ref2 = RevBayesCore::RbHelpReference();
      ref2.setCitation("Algorithmic improvements to species delimitation and phylogeny estimation under the multispecies coalescent. Graham Jones.  Journal of Mathematical Biology, 2016.");
      ref2.setDoi("https://doi.org/10.1007/s00285-016-1034-0");
      ref2.setUrl("http://link.springer.com/article/10.1007/s00285-016-1034-0");

      references.push_back(ref2);

      return references;
  }


  /**
   * Get the names of similar and suggested other functions
   */
  std::vector<std::string> Move_SpeciesSubtreeScale::getHelpSeeAlso(void) const
  {
      // create an entry for each suggested function
      std::vector<std::string> see_also;
      see_also.push_back( "mvSpeciesNodeTimeSlideUniform" );
      see_also.push_back( "mvSpeciesSubtreeScaleBeta" );
      see_also.push_back( "mvSpeciesNarrow" );
      see_also.push_back( "mvSpeciesTreeScale" );

      return see_also;
  }


  /**
   * Get the title of this help entry
   */
  std::string Move_SpeciesSubtreeScale::getHelpTitle(void) const
  {
      // create a title variable
      std::string title = "Subtree scale move on species tree and gene trees for multispecies coalescent models. ";

      return title;
  }
