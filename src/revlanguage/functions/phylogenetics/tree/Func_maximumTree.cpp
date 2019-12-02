#include "Func_maximumTree.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlTimeTree.h"
#include "ModelVector.h"
#include "MaximumTreeFunction.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_maximumTree::Func_maximumTree( void ) : TypedFunction<TimeTree>( )
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_maximumTree* Func_maximumTree::clone( void ) const
{

    return new Func_maximumTree( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::Tree >* Func_maximumTree::createFunction( void ) const
{

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::Tree > >* gTrees = static_cast<const ModelVector< TimeTree > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::MaximumTreeFunction* f = new RevBayesCore::MaximumTreeFunction( gTrees );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_maximumTree::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {

        argumentRules.push_back( new ArgumentRule( "geneTrees", ModelVector< TimeTree >::getClassTypeSpec() , "The vector of trees from which to pick the maximum.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }

    return argumentRules;
}


const std::string& Func_maximumTree::getClassType(void)
{

    static std::string rev_type = "Func_maximumTree";

	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_maximumTree::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_maximumTree::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "maximumTree";

    return f_name;
}


const TypeSpec& Func_maximumTree::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_maximumTree::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Bastien Boussau" );

    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_maximumTree::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Builds the maximum species tree given several ultrametric gene trees.";

    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_maximumTree::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details =  "The maximum species tree is a consistent estimate of the species tree under the multispecies coalescent model, if the gene trees are correct and the effective population size constant along the species tree." ;

    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_maximumTree::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "# We simulate a species tree and gene trees and reconstruct a species tree using maximum tree:\n";
    example += "# Let’s simulate a species tree with 10 taxa, 5 gene trees, 1 alleles per species:\n";
    example += "n_species <- 10\n";
    example += "n_genes <- 5\n";
    example += "n_alleles <- 2\n";
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
    example += "        taxons[g][(i-1)*n_alleles+j] <- taxon(taxonName=\"Species_\"+i, speciesName=\"Species_\"+i)\n";
    example += "    }\n";
    example += "  }\n";
    example += "  geneTrees[g] ~ dnMultiSpeciesCoalescent(speciesTree=spTree, Ne=popSize, taxa=taxons[g])\n";
    example += "}\n";
    example += "# Let's compute the maximum tree: \n";
    example += "recTree <- maximumTree(geneTrees)\n";
    example += "print(recTree)\n";


    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_maximumTree::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    RevBayesCore::RbHelpReference ref = RevBayesCore::RbHelpReference();
    ref.setCitation("High-resolution species trees without concatenation. Scott V. Edwards, Liang Liu, and Dennis K. Pearl . PNAS April 3, 2007 vol. 104 no. 14 .");
    ref.setDoi("");
    ref.setUrl("http://www.pnas.org/content/104/14/5936.full");

    references.push_back(ref);

    RevBayesCore::RbHelpReference ref2 = RevBayesCore::RbHelpReference();
    ref2.setCitation("Maximum tree: a consistent estimator of the species tree. Liu L, Yu L, Pearl DK.  Journal of Mathematical Biology, 2010. Jan;60(1):95-106.");
    ref2.setDoi("https://doi.org/10.1007/s00285-009-0260-0");
    ref2.setUrl("https://link.springer.com/article/10.1007%2Fs00285-009-0260-0");

    references.push_back(ref2);



    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_maximumTree::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;


    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_maximumTree::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Maximum tree function to build a species tree.";

    return title;
}
