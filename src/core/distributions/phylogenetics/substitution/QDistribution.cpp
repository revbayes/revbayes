#include "QDistribution.h"

#include <cmath>
#include "RateMatrix_MPQ.h"
#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;



QDistribution::QDistribution (const TypedDagNode<RbVector<double> >* a, double lrr, double lrnr) : 
    TypedDistribution<RateGenerator>(new RateMatrix_MPQ()), alpha(a), log_rho_reversible(lrr), log_rho_non_reversible(lrnr) {

    addParameter(alpha);
    
    if (alpha->getValue().size() != 4)
        {
        throw RbException("You need to provide 4 elements for the prior.");
        }
    redrawValue();
}

QDistribution* QDistribution::clone(void) const {

    return new QDistribution( *this );
}

double QDistribution::computeLnProbability(void) {

    RateMatrix_MPQ& myRateMatrix = static_cast<RateMatrix_MPQ&>(*this->value);
    
    /* Read reversibility from the rate matrix itself.  This class used to keep
       its own isReversible flag, which was never assigned anywhere, so the
       polytope volume below was being chosen by an uninitialized bool. */
    bool is_reversible = myRateMatrix.getIsReversible();

    // get the current value for alpha (prior on stationary frequencies)
    const std::vector<double>& currentAlpha = alpha->getValue();

    // calcualte the prior density for the base frequencies
    const std::vector<double>& bf = myRateMatrix.getStationaryFrequencies();
    double lnProb = RbStatistics::Dirichlet::lnPdf(currentAlpha, bf);
    
    // the prior on the model itself
    if (is_reversible == true)
        lnProb += log_rho_reversible;
    else
        lnProb += log_rho_non_reversible;

    /* Given pi, the weights w_ij = pi_i q_ij are uniform on the polytope of the
       current model, so the density is the reciprocal of that polytope's volume
       measured in the same free coordinates the reversible-jump Jacobian is taken
       in: (w_AC, w_AG, w_AT, w_CG, w_CT) for the time-reversible model, a 5-simplex
       of volume 1/(5! 2^5) = 1/3840, and (w_AC, w_AG, w_AT, w_CA, w_CG, w_CT, w_GA,
       w_GC) for the non-reversible model, the circulation polytope of K4 with volume
       1/(48 8!) = 1/1935360.  Both are constants.  There is no dependence on the
       stationary frequencies and none on the rates: the pairing of a flat prior in
       w with a Jacobian taken in w is what makes the acceptance probability well
       defined, and a mismatched pairing is easy to introduce and hard to notice. */
    if (is_reversible == true)
        lnProb += std::log(3840.0);
    else
        lnProb += std::log(1935360.0);

    return lnProb;
}

/* Set the prior odds log(rho_R / rho_N) while keeping rho_R + rho_N = 1.

   Only the difference of the two matters to the MCMC, so the normalization is
   for reporting, but it has to be done in a form that survives a large tilt:
   the whole point of tuning is that delta may have to reach into the tens or
   hundreds before a strongly supported model stops monopolizing the chain, and
   exp(delta) overflows long before log1p of it does. */
void QDistribution::setLnPriorOdds(double delta) {

    if (delta > 0.0)
        {
        log_rho_reversible     = -log1p(exp(-delta));
        log_rho_non_reversible = -delta - log1p(exp(-delta));
        }
    else
        {
        log_rho_reversible     = delta - log1p(exp(delta));
        log_rho_non_reversible = -log1p(exp(delta));
        }
}

void QDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, Boolean &rv) const {

    if ( n == "isReversible" )
        {
        RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);
        rv = my_rate_matrix.getIsReversible();
        }
    else
        {
        throw RbException("The Q-Distribution does not have a member method called '" + n + "'.");
        }
}

void QDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, Simplex &rv) const {

    if (n == "getPi")
        {
        RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);
        rv = my_rate_matrix.getStationaryFrequencies();
        }
    else
        {
        throw RbException("The Q-Distribution does not have a member method called '" + n + "'.");
        }
}

void QDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const {

    if (n == "getRates")
        {
        RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);
        rv = my_rate_matrix.getRates();
        }
    else
        {
        throw RbException("The Q-Distribution does not have a member method called '" + n + "'.");
        }
}


/* The prior on the model indicator, exposed so that it can be monitored.

   This has to be logged whenever the prior is tuned, because the sampled model
   frequencies mean nothing on their own: the Bayes factor is recovered from them
   together with the prior odds that produced them,

       BF_NR = (p_N / p_R) * exp(lnPriorOdds)

   and the posterior probability under equal priors is BF_NR / (1 + BF_NR). */
void QDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, double &rv) const {

    if ( n == "lnPriorOdds" )
        {
        rv = log_rho_reversible - log_rho_non_reversible;
        }
    else if ( n == "lnRhoReversible" )
        {
        rv = log_rho_reversible;
        }
    else if ( n == "lnRhoNonReversible" )
        {
        rv = log_rho_non_reversible;
        }
    else
        {
        throw RbException("The Q-Distribution does not have a member method called '" + n + "'.");
        }
}


void QDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == alpha)
        {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
        }
}


//void QDistribution::keepSpecialization(DagNode *toucher);
//void QDistribution::restoreSpecialization(DagNode *toucher);
//void QDistribution::touchSpecialization(DagNode *toucher, bool touchAll);

void QDistribution::redrawValue(void) {
    
    // implement how to draw a new value from the prior
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);

    const std::vector<double>& curr_alpha = alpha->getValue();
    
    double u = log( rng->uniform01() );
    
    if ( u < log_rho_reversible )
        {
        my_rate_matrix.initializeTimeReversibleModel(curr_alpha, rng);
        }
    else
        {
        my_rate_matrix.initializeNonReversibleModel(curr_alpha, rng);
        }
    
    my_rate_matrix.update();
}
