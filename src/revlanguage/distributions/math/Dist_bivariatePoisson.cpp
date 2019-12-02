#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Dist_bivariatePoisson.h"
#include "BivariatePoissonDistribution.h"
#include "Natural.h"
#include "ModelVector.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "StochasticNode.h"

using namespace RevLanguage;


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_bivariatePoisson::Dist_bivariatePoisson() : TypedDistribution< ModelVector<Natural> >()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
Dist_bivariatePoisson* Dist_bivariatePoisson::clone( void ) const
{

    return new Dist_bivariatePoisson(*this);
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
RevBayesCore::BivariatePoissonDistribution* Dist_bivariatePoisson::createDistribution( void ) const
{
    // get the parameters
  RevBayesCore::TypedDagNode<double>* theta1x = static_cast<const RealPos &>( theta1->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode<double>* theta2x = static_cast<const RealPos &>( theta2->getRevObject() ).getDagNode();
  RevBayesCore::TypedDagNode<double>* theta0x = static_cast<const RealPos &>( theta0->getRevObject() ).getDagNode();
  RevBayesCore::BivariatePoissonDistribution* d        = new RevBayesCore::BivariatePoissonDistribution( theta1x, theta2x, theta0x );

    return d;
}



/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_bivariatePoisson::getClassType(void)
{

    static std::string rev_type = "Dist_bivariatePoisson";

	return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_bivariatePoisson::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Natural>::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Dist_bivariatePoisson::getDistributionFunctionAliases( void ) const
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
std::string Dist_bivariatePoisson::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "bivariatePoisson";

    return d_name;
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Dist_bivariatePoisson::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Alexander Zarebski" );
    authors.push_back( "Sebastian Hoehna" );

    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Dist_bivariatePoisson::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A Bivariate Poisson distribution defines probabilities for pairs of natural numbers.";

    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Dist_bivariatePoisson::getHelpDetails(void) const
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
std::string Dist_bivariatePoisson::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";

    example += "th1 ~ dnUniform(0.0,10.0)\n";
    example += "th2 ~ dnUniform(0.0,10.0)\n";
    example += "th0 ~ dnUniform(0.0,10.0)\n";
    example += " ~ dnBivariatePoisson(th1, th2, th0)\n";
    example += "x.clamp([3, 3, 3])\n";
    example += "oves[1] = mvSlide(th1, delta=0.01, weight=1.0)\n";
    example += "moves[2] = mvSlide(th2, delta=0.01, weight=1.0)\n";
    example += "oves[3] = mvSlide(th0, delta=0.01, weight=1.0)\n";
    example += "monitors[1] = mnScreen(printgen=1000,  th0)\n";
    example += "ymodel = model(th1)\n";
    example += "mymcmc = mcmc(mymodel, monitors, moves)\n";
    example += "ymcmc.burnin(generations=20000,tuningInterval=100)\n";
    example += "mymcmc.run(generations=200000)\n";

    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Dist_bivariatePoisson::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;

    RevBayesCore::RbHelpReference ref1;
    ref1.setCitation("Karlis D, Ntzoufras J (2003). Bayesian and Non-Bayesian Analysis of Soccer Data using Bivariate Poisson Regression Models.");
    references.push_back( ref1 );

    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Dist_bivariatePoisson::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "dnPoisson" );


    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Dist_bivariatePoisson::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Bivariate Poisson Distribution";

    return title;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the Poisson distribution are:
 * (1) the rate lambda which must be a positive real between 0 and 1 (= a probability).
 *
 * \return The member rules.
 */
const MemberRules& Dist_bivariatePoisson::getParameterRules(void) const
{

    static MemberRules distBivPoisMemberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
      distBivPoisMemberRules.push_back( new ArgumentRule( "theta1", RealPos::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ));
      distBivPoisMemberRules.push_back( new ArgumentRule( "theta2", RealPos::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ));
      distBivPoisMemberRules.push_back( new ArgumentRule( "theta0", RealPos::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ));

        rules_set = true;
    }

    return distBivPoisMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_bivariatePoisson::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
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
void Dist_bivariatePoisson::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

  if ( name == "theta1" ) {
    theta1 = var;
  } else if (name == "theta2") {
    theta2 = var;
  } else if (name == "theta0") {
    theta0 = var;
  }
    else
    {
      TypedDistribution< ModelVector<Natural> >::setConstParameter(name, var);
    }

}
