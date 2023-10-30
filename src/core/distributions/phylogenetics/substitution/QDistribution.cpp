#include "QDistribution.h"
#include "RateMatrix_MPQ.h"
#include "DistributionDirichlet.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;


QDistribution::QDistribution (const TypedDagNode<RbVector<double> >* a,
                              const TypedDagNode<double>* r) : TypedDistribution<RateGenerator>( new RateMatrix_MPQ() ),
alpha( a ),
rho( r )
{
    addParameter( alpha );
    addParameter( rho );
    
    if ( alpha->getValue().size() != 4 )
    {
        throw RbException("You need to provide 4 elements for the prior.");
    }
    
    redrawValue();
    
}



QDistribution* QDistribution::clone(void) const
{
    return new QDistribution( *this );
}


double QDistribution::computeLnProbability(void)
{
    RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);
    
    // get the current value for alpha (prior on stationary frequencies)
    const std::vector<double>& curr_alpha = alpha->getValue();
    double model_prob = rho->getValue();
    
    // compute the prior on the stationary frequencies
    std::vector<mpq_class>& pi = my_rate_matrix.getPi();
    std::vector<double> bf(4);
    for (int i=0; i<4; i++)
        bf[i] = pi[i].get_d();
    double ln_prob = RbStatistics::Dirichlet::lnPdf(curr_alpha, bf);
    
    // compute prior of the rate matrix
    ln_prob += 3.0 * (log(bf[0]) + log(bf[1]) + log(bf[2]) + log(bf[3]));
    if (my_rate_matrix.getIsReversible() == true)
        {
        ln_prob += log(120.0);
        ln_prob += log(4.0) + 0.5 * log(3.0);
        ln_prob += log( model_prob );
        }
    else
        {
        ln_prob += log(40320.0);
        ln_prob -= 0.5 * log(6.0);
        ln_prob += log( 1.0 - model_prob );
        }

    
    return ln_prob;
}


void QDistribution::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, Boolean &rv) const
{

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


void QDistribution::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< double >* >( newP );
    }
    
}


//void QDistribution::keepSpecialization(DagNode *toucher);
//void QDistribution::restoreSpecialization(DagNode *toucher);
//void QDistribution::touchSpecialization(DagNode *toucher, bool touchAll);

void QDistribution::redrawValue(void)
{
    
    // implement how to draw a new value from the prior
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    RateMatrix_MPQ& my_rate_matrix = static_cast<RateMatrix_MPQ&>(*this->value);
    
    double model_prob = rho->getValue();
    const std::vector<double>& curr_alpha = alpha->getValue();
    
    double u = rng->uniform01();
    if ( u < model_prob )
    {
        my_rate_matrix.initializeTimeReversibleModel(curr_alpha, rng);
    }
    else
    {
        my_rate_matrix.initializeNonReversibleModel(curr_alpha, rng);
    }
    
    my_rate_matrix.update();
}
