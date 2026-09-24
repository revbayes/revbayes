#include <math.h>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "BetaDistribution.h"
#include "Dist_Q.h"
#include "RbException.h"
#include "Real.h"
#include "RlString.h"
#include "RlSimplex.h"
#include "Probability.h"
#include "RlContinuousStochasticNode.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlDistributionMemberFunction.h"
#include "RlStochasticNode.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { class ContinuousDistribution; }

using namespace RevLanguage;

Dist_Q::Dist_Q(void) : TypedDistribution<RateGenerator>() {
    
    setGuiDistributionName("Q");
    setGuiDistributionToolTip("Q distribution for random variables on rate matrices");
}

Dist_Q::~Dist_Q(void) {
    
}

Dist_Q* Dist_Q::clone(void) const {
    
    return new Dist_Q(*this);
}

RevBayesCore::QDistribution* Dist_Q::createDistribution(void) const {
    
    // get the parameters
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* a = static_cast<const ModelVector<RealPos> &>(alpha->getRevObject()).getDagNode();
        
    double logRhoReversible    = 0.0;
    double logRhoNonReversible = 0.0;

    double delta = static_cast<const Real &>( lnPriorOdds->getRevObject() ).getDagNode()->getValue();

    if (delta != 0.0)
        {
        /* The prior odds given directly, as log(rho_R / rho_N).

           This route exists because rho cannot express a large tilt. rho is
           log(rho_R), and log(rho_N) has to be recovered from it as
           log1p(-exp(rho)); but for a tilt beyond about 36, rho_R differs from one
           by less than a double can represent, exp(rho) is exactly 1, and the
           recovery is log(0). Working with the odds instead keeps both weights
           accurate at any magnitude, which matters when a fixed prior is being set
           deliberately far from even. */
        if (delta > 0.0)
            {
            logRhoReversible    = -log1p( exp(-delta) );
            logRhoNonReversible = -delta - log1p( exp(-delta) );
            }
        else
            {
            logRhoReversible    = delta - log1p( exp(delta) );
            logRhoNonReversible = -log1p( exp(delta) );
            }
        }
    else
        {
        /* rho is the LOG of the prior probability of the time-reversible model, so it
           must be negative. Passing a probability rather than its log, say rho=0.5,
           silently produced log(1 - exp(0.5)) = log of a negative number and seeded
           the whole analysis with a NaN model prior, which then made every acceptance
           ratio NaN. Catch it here instead. */
        double tmp = static_cast<const Real &>( rho->getRevObject() ).getDagNode()->getValue();
        if ( tmp >= 0.0 )
            {
            throw RbException("The 'rho' argument of dnQ is the LOG of the prior probability of the time-reversible model, so it must be negative. Use rho=ln(0.5) rather than rho=0.5.");
            }
        if ( tmp < -36.0 )
            {
            throw RbException("The 'rho' argument of dnQ cannot express a prior tilt beyond about 36 log units, because log(rho_N) has to be recovered from it and the recovery underflows. Use the 'lnPriorOdds' argument instead.");
            }
        logRhoReversible    = tmp;
        logRhoNonReversible = log1p( -exp(logRhoReversible) );
        }
        
    /* Confining the chain to one model is how the two are measured separately, which
       is the only way to compare their likelihoods with none of the reversible-jump
       machinery involved. It also switches the jump moves off, because
       MPQRateMatrixProposal asks the distribution rather than keeping its own flag. */
    const std::string& fm = static_cast<const RlString &>(fixedModel->getRevObject()).getValue();
    RevBayesCore::QDistribution::FIXED_MODEL fixed = RevBayesCore::QDistribution::FREE;
    if ( fm == "reversible" )
        fixed = RevBayesCore::QDistribution::REVERSIBLE_ONLY;
    else if ( fm == "nonreversible" )
        fixed = RevBayesCore::QDistribution::NON_REVERSIBLE_ONLY;
    else if ( fm != "free" )
        throw RbException("The 'fixedModel' argument of dnQ must be \"free\", \"reversible\" or \"nonreversible\".");

    RevBayesCore::QDistribution* d = new RevBayesCore::QDistribution(a, logRhoReversible, logRhoNonReversible, fixed);
    return d;
}

RateMatrix* Dist_Q::createRandomVariable(void) const {
    
    RevBayesCore::TypedDistribution<RevBayesCore::RateGenerator>* d = createDistribution();
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* rv  = new StochasticNode("", d, this->clone() );
    return new RateMatrix(rv);
}

/* Get Rev type of object */
const std::string& Dist_Q::getClassType(void) {
    
    static std::string rev_type = "Dist_Q";
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Dist_Q::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<Probability>::getClassTypeSpec() ) );
    return rev_type_spec;
}

/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * \return Rev name of constructor function.
 */
std::string Dist_Q::getDistributionFunctionName(void) const {

    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "Q";
    return d_name;
}

MethodTable Dist_Q::getDistributionMethods(void) const {
    
    MethodTable methods = TypedDistribution< RateGenerator >::getDistributionMethods();
    
    // member functions
    ArgumentRules* get_reversibility_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_Q, RlBoolean >( "isReversible", this->variable, get_reversibility_arg_rules, true ) );

    ArgumentRules* get_pi_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_Q, Simplex >( "getPi", this->variable, get_pi_arg_rules, true ) );

    ArgumentRules* get_rates_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_Q, ModelVector<RealPos> >( "getRates", this->variable, get_rates_arg_rules, true ) );

    /* Where in the polyhedron the non-reversible model sits: the three coordinates
       are one half exactly at time reversibility. Monitor these to find out whether
       the polyhedron draw of mvMPQRateMatrix should be concentrated toward the center
       or spread toward the facets, which is what polyhedronAlpha controls. */
    ArgumentRules* get_u_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_Q, ModelVector<RealPos> >( "getU", this->variable, get_u_arg_rules, true ) );

    /* The prior on the model indicator. These matter when the prior is being
       tuned by mvMPQRateMatrix(tuneModelPrior=TRUE), because the sampled model
       frequencies can only be turned back into a Bayes factor if the prior odds
       that produced them are known:

           BF_NR = (p_N / p_R) * exp(lnPriorOdds)

       Monitor lnPriorOdds alongside isReversible and the tuned value is on record
       with the samples it produced. */
    ArgumentRules* get_ln_prior_odds_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_Q, Real >( "lnPriorOdds", this->variable, get_ln_prior_odds_arg_rules, true ) );

    ArgumentRules* get_ln_rho_rev_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_Q, Real >( "lnRhoReversible", this->variable, get_ln_rho_rev_arg_rules, true ) );

    ArgumentRules* get_ln_rho_nonrev_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_Q, Real >( "lnRhoNonReversible", this->variable, get_ln_rho_nonrev_arg_rules, true ) );

    return methods;
}

/** Return member rules (no members) */
const MemberRules& Dist_Q::getParameterRules(void) const {
    
    static MemberRules distQMemberRules;
    static bool rulesSet = false;
    
    if (!rulesSet)
        {
        distQMemberRules.push_back( new ArgumentRule( "alpha",  ModelVector<RealPos>::getClassTypeSpec(), "The alpha shape parameter.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        distQMemberRules.push_back( new ArgumentRule( "rho" , Real::getClassTypeSpec(), "The log probability of the reversible model. Must be negative, and cannot express a tilt beyond about 36 log units; use lnPriorOdds for anything larger.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real( -0.6931471805599453 ) ) );
        distQMemberRules.push_back( new ArgumentRule( "fixedModel" , RlString::getClassTypeSpec(), "Confine the chain to one model: \"free\" (the default) lets the reversible jump operate, while \"reversible\" or \"nonreversible\" starts the chain in that model and suppresses the jump moves entirely. Use this to measure the two models separately.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString( "free" ) ) );
        distQMemberRules.push_back( new ArgumentRule( "lnPriorOdds" , Real::getClassTypeSpec(), "The prior odds on the model indicator, given directly as log(rho_R/rho_N). When this is nonzero it is used and 'rho' is ignored. Unlike 'rho' it stays accurate at any magnitude, so it is the argument to use when setting a deliberately extreme prior.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real( 0.0 ) ) );
        rulesSet = true;
        }
    
    return distQMemberRules;
}

const TypeSpec& Dist_Q::getTypeSpec(void) const {
    
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}

/** Print value for user */
void Dist_Q::printValue(std::ostream& o) const {
    
    o << "Q(alpha=";
    if ( alpha != NULL )
        {
        o << alpha->getName();
        }
    else
        {
        o << "?";
        }
    o << ", prob=";
    if ( rho != NULL )
        {
        o << rho->getName();
        }
    else
        {
        o << "?";
        }
    o << ")";
}

/** Set a member variable */
void Dist_Q::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
        
    if (name == "alpha")
        {
        alpha = var;
        }
    else if (name == "rho")
        {
        rho = var;
        }
    else if (name == "lnPriorOdds")
        {
        lnPriorOdds = var;
        }
    else if (name == "fixedModel")
        {
        fixedModel = var;
        }
    else
        {
        TypedDistribution<RateGenerator>::setConstParameter(name, var);
        }
}
