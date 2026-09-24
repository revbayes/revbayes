#include "QDistribution.h"

#include <cmath>
#include <iostream>
#include "RateMatrix_MPQ.h"
#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;



QDistribution::QDistribution (const TypedDagNode<RbVector<double> >* a, double lrr, double lrnr, FIXED_MODEL fm) : 
    TypedDistribution<RateGenerator>(new RateMatrix_MPQ()), alpha(a), log_rho_reversible(lrr), log_rho_non_reversible(lrnr), fixed_model(fm) {

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

    if (is_reversible == true)
        lnProb += std::log(3840.0);
    else
        lnProb += std::log(1935360.0);

#   ifdef MPQ_DEBUG_PRIOR
        {
        static long n_rev = 0;
        static long n_non = 0;
        long& n = ( is_reversible ? n_rev : n_non );
        n++;
        if ( n <= 3 )
            {
            std::cerr << "QDistribution prior [" << (is_reversible ? "reversible" : "non-reversible") << "]"
                      << "  Dirichlet = " << RbStatistics::Dirichlet::lnPdf(currentAlpha, bf)
                      << "  lnRho = "     << (is_reversible ? log_rho_reversible : log_rho_non_reversible)
                      << "  lnInvVolume = " << (is_reversible ? std::log(3840.0) : std::log(1935360.0))
                      << "  total = "     << lnProb << std::endl;
            std::cerr << "    alpha =";
            for (size_t i=0; i<currentAlpha.size(); i++)
                std::cerr << " " << currentAlpha[i];
            std::cerr << "    pi =";
            for (size_t i=0; i<bf.size(); i++)
                std::cerr << " " << bf[i];
            std::cerr << std::endl;
            }
        }
#   endif

    return lnProb;
}

/* Set the prior odds log(rho_R / rho_N) while keeping rho_R + rho_N = 1 */
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
    else if (n == "getU")
        {
        // where in the polyhedron the current non-reversible model sits.
        RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);
        mpq_class u1, u2, u3;
        rv.clear();
        if ( my_rate_matrix.getIsReversible() == true )
            {
            rv.push_back( 0.5 );
            rv.push_back( 0.5 );
            rv.push_back( 0.5 );
            }
        else if ( my_rate_matrix.recoverU(u1, u2, u3) == true )
            {
            rv.push_back( u1.get_d() );
            rv.push_back( u2.get_d() );
            rv.push_back( u3.get_d() );
            }
        else
            {
            throw RbException("Could not recover the point in the polyhedron from the rate matrix.");
            }
        }
    else
        {
        throw RbException("The Q-Distribution does not have a member method called '" + n + "'.");
        }
}


/* The prior on the model indicator, exposed so that it can be monitored.

   This has to be logged whenever the prior is tuned, because the sampled model
   frequencies mean nothing on their own. The Bayes factor is recovered from them
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

void QDistribution::redrawValue(void) {
    
    // implement how to draw a new value from the prior
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);

    const std::vector<double>& curr_alpha = alpha->getValue();
    
    /* When the distribution is confined to one model, start there; otherwise draw
       the model from its prior. */
    bool start_reversible;
    if ( fixed_model == REVERSIBLE_ONLY )
        start_reversible = true;
    else if ( fixed_model == NON_REVERSIBLE_ONLY )
        start_reversible = false;
    else
        start_reversible = ( log( rng->uniform01() ) < log_rho_reversible );

    if ( start_reversible == true )
        {
        my_rate_matrix.initializeTimeReversibleModel(curr_alpha, rng);
        }
    else
        {
        my_rate_matrix.initializeNonReversibleModel(curr_alpha, rng);
        }
    
    my_rate_matrix.update();
}
