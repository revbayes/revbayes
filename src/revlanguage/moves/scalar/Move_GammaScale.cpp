#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "ContinuousStochasticNode.h"
#include "MetropolisHastingsMove.h"
#include "Move_GammaScale.h"
#include "Probability.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "GammaScaleProposal.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_GammaScale::Move_GammaScale() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_GammaScale* Move_GammaScale::clone(void) const
{
    
    return new Move_GammaScale(*this);
}


/**
 * Create a new internal move object.
 *
 * This function simply dynamically allocates a new internal move object that is
 * associated with the variable (DAG-node). The internal move object is created by calling its
 * constructor and passing the move-parameters (the variable and other parameters) as arguments of the
 * constructor. The move constructor takes care of the proper hook-ups.
 *
 */
void Move_GammaScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    RevBayesCore::Proposal *p = NULL;
    
    // now allocate a new sliding move
    double d = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double r = static_cast<const Probability &>( tuneTarget->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<double>* tmp = static_cast<const RealPos &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<double> *n = dynamic_cast<RevBayesCore::StochasticNode<double> *>( tmp );
    p = new RevBayesCore::GammaScaleProposal(n, d, r);
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    value = new RevBayesCore::MetropolisHastingsMove(p, w, t);
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_GammaScale::getClassType(void)
{
    
    static std::string rev_type = "Move_GammaScale";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_GammaScale::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_GammaScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "GammaScale";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the gamma scale move are:
 * (1) the variable must be a positive real.
 * (2) the tuning parameter lambda that defines the size of the proposal must be a positive real
 * (3) a flag whether auto-tuning should be used.
 *
 * \return The member rules.
 */
const MemberRules& Move_GammaScale::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"     , RealPos::getClassTypeSpec()  , "The variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()  , "The strength of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(), "Should we tune lambda during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_GammaScale::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



void Move_GammaScale::printValue(std::ostream &o) const
{
    
    o << "GammaScale(";
    if (x != NULL)
    {
        o << x->getName();
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
void Move_GammaScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "tune" )
    {
        tune = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}

/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Move_GammaScale::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Jeremy M. Brown" );
    return authors;
}

/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Move_GammaScale::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "mvScale" );
    return see_also;
}

/**
 * Get the (brief) description for this move
 */
std::string Move_GammaScale::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "A move to scale a single continuous value by multiplying by a value drawn from a Gamma(lambda,1) distribution. Lambda is the tuning parameter that controls the size of the proposals.";
    return description;
}

/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Move_GammaScale::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    
    example += "# Here is a simple example for conducting MCMC on the mean and sd of a Normal distribution.\n";
    example += "\n";
    example += "# Uniform(0,1) priors on the mean and sd\n";
    example += "mean ~ dnUnif(0,1)\n";
    example += "sd ~ dnUnif(0,1)\n";
    example += "\n";
    example += "# Dummy data (will not actually be analyzed)\n";
    example += "data <- v(0.4,0.5,0.6)\n";
    example += "\n";
    example += "# Clamping data\n";
    example += "for (i in 1:data.size()){ outcomes[i] ~ dnNorm(mean,sd); outcomes[i].clamp(data[i]) }\n";
    example += "\n";
    example += "# Initializing move and monitor counters\n";
    example += "mvi = 1\n";
    example += "mni = 1\n";
    example += "\n";
    example += "# Adding Gamma scale moves for the mean and sd (THIS MOVE IS HERE)\n";
    example += "moves[mvi++] = mvGammaScale(mean)\n";
    example += "moves[mvi++] = mvGammaScale(sd)\n";
    example += "\n";
    example += "# Instantiating the model\n";
    example += "mymodel = model(outcomes)\n";
    example += "\n";
    example += "# Adding screen monitor for the mean\n";
    example += "monitors[mni++] = mnScreen(mean, printgen=1000)\n";
    example += "\n";
    example += "# Creating MCMC object\n";
    example += "mymcmc = mcmc(mymodel, moves, monitors)\n";
    example += "\n";
    example += "# Running MCMC under the prior\n";
    example += "mymcmc.run(30000,underPrior=TRUE);\n";
    
    return example;
}


