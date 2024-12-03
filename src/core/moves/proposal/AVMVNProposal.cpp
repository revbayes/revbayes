#include "AVMVNProposal.h"

#include <cmath>
#include <algorithm>

#include "CholeskyDecomposition.h"
#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RbException.h"
#include "ContinuousStochasticNode.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Simplex.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class RandomNumberGenerator; }

using namespace RevBayesCore;

/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
AVMVNProposal::AVMVNProposal( double s, double e, size_t n0, size_t c0, size_t m ) : Proposal(),
    noTransformScalarVariables(  ),
    logTransformScalarVariables(  ),
    logitTransformScalarVariables(  ),
    noTransformVectorVariables(  ),
    logConstrainedSumTransformVectorVariables(  ),
    lnHastingsratio( 0.0 ),
    sigma( s ),
    epsilon( e ),
    waitBeforeLearning( n0 ),
    waitBeforeUsing ( c0 ),
    maxUpdates ( m ),
    nTried ( 0 ),
    updates ( 0 ),
    dim ( 0 ),
    C_emp(  ),
    AVMVN_cholesky_L(  ),
    x_bar(  ),
    storedValues(  ),
    storedValuesUntransformed(  ),
    proposedValues(  )
{
    if (waitBeforeUsing < waitBeforeLearning)
    {
        throw RbException("Cannot delay learning empirical covariance matrix longer than using empirical covariance matrix");
    }
}


/**
 * Add a variable.
 */
void AVMVNProposal::addUntransformedScalar( StochasticNode<double> *v )
{
    // add it to the nodes vector
    addNode( v );

    noTransformScalarVariables.push_back( v );
}

void AVMVNProposal::addLogScalar( StochasticNode<double> *v )
{
    // add it to the nodes vector
    addNode( v );

    logTransformScalarVariables.push_back( v );
}

void AVMVNProposal::addLogitScalar( ContinuousStochasticNode *v )
{
    // add it to the nodes vector
    addNode( v );

    logitTransformScalarVariables.push_back( v );
}

/**
 * Add a variable.
 */

void AVMVNProposal::addUntransformedVector( StochasticNode<RbVector<double> > *v)
{
    // add it to the nodes vector
    addNode( v );

    noTransformVectorVariables.push_back( v );
}

void AVMVNProposal::addLogConstrainedSumVector( StochasticNode<Simplex> *v)
{
    // add it to the nodes vector
    addNode( v );

    logConstrainedSumTransformVectorVariables.push_back( v );

}

/**
 * Calculates hastings ratio of the move, also sets all variables to their new values.
 */
void AVMVNProposal::calculateHastingsRatio( std::vector<double> x_prime, std::vector<double> x )
{
    // Calculate hastings ratio and give variables their new values
    size_t index = 0;
    for (size_t i = 0; i < noTransformScalarVariables.size(); ++i)
    {
        noTransformScalarVariables[i]->getValue() = x_prime[index];
        ++index;
    }

    for (size_t i = 0; i < logTransformScalarVariables.size(); ++i)
    {
        logTransformScalarVariables[i]->getValue() = exp(x_prime[index]);
        lnHastingsratio += -x[index] + x_prime[index];
        ++index;
    }

    for (size_t i = 0; i < logitTransformScalarVariables.size(); ++i)
    {
        double p_old = 1 / (1 + exp(-x[index]));
        double p_proposed = 1 / (1 + exp(-x_prime[index]));

        double lb = logitTransformScalarVariables[i]->getMin();
        double ub = logitTransformScalarVariables[i]->getMax();

        logitTransformScalarVariables[i]->getValue() = p_proposed * (ub - lb) + lb;
        lnHastingsratio += ( -log(1.0 - p_old) - log(p_old) ) - ( -log(1.0 - p_proposed) - log(p_proposed) );
        ++index;
    }

    for (size_t i = 0; i < noTransformVectorVariables.size(); ++i)
    {
        std::vector<double> tmp;

        for (size_t j=0; j < noTransformVectorVariables[i]->getValue().size(); ++j)
        {
            tmp.push_back( x_prime[index] );
            ++index;
        }
        noTransformVectorVariables[i]->getValue() = tmp;
    }

    for (size_t i = 0; i < logConstrainedSumTransformVectorVariables.size(); ++i)
    {
        std::vector<double> old_sum_constrained = logConstrainedSumTransformVectorVariables[i]->getValue();
        std::vector<double> tmp;
        double sum = 0.0;

        for (size_t j=0; j < old_sum_constrained.size(); ++j)
        {
            tmp.push_back( exp(x_prime[index]) );
            sum += tmp[j];
            ++index;
        }
        for (size_t j=0; j < old_sum_constrained.size(); ++j)
        {
            tmp[j] /= sum;
            lnHastingsratio += -log(old_sum_constrained[j]) + log(tmp[j]);
        }

        logConstrainedSumTransformVectorVariables[i]->getValue() = tmp;

    }
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
void AVMVNProposal::cleanProposal( void )
{

}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
AVMVNProposal* AVMVNProposal::clone( void ) const
{

    return new AVMVNProposal( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
const std::string& AVMVNProposal::getProposalName( void ) const
{
    static std::string name = "AVMVN";

    return name;
}


double AVMVNProposal::getProposalTuningParameter( void ) const
{
    return sigma;
}


/**
 * Obtains the values of all variables in the AVMVN proposal on the transformed scale.
 */
void AVMVNProposal::getAVMVNMemberVariableValues( std::vector<double> *x )
{
    for (size_t i = 0; i < noTransformScalarVariables.size(); ++i)
    {
        x->push_back( noTransformScalarVariables[i]->getValue() );
    }

    for (size_t i = 0; i < logTransformScalarVariables.size(); ++i)
    {
        x->push_back( log(logTransformScalarVariables[i]->getValue()) );
    }

    for (size_t i = 0; i < logitTransformScalarVariables.size(); ++i)
    {
        double y = logitTransformScalarVariables[i]->getValue();
        double lb = logitTransformScalarVariables[i]->getMin();
        double ub = logitTransformScalarVariables[i]->getMax();
        double p = (y - lb) / (ub - lb);
        x->push_back( log(p / (1 - p)) );

    }

    for (size_t i = 0; i < noTransformVectorVariables.size(); ++i)
    {

        std::vector<double> tmp = noTransformVectorVariables[i]->getValue();

        for (size_t j=0; j < tmp.size(); ++j)
        {
            x->push_back( tmp[j] );
        }
    }

    for (size_t i = 0; i < logConstrainedSumTransformVectorVariables.size(); ++i)
    {
        std::vector<double> tmp = logConstrainedSumTransformVectorVariables[i]->getValue();

        for (size_t j=0; j < tmp.size(); ++j)
        {
            x->push_back( log(tmp[j]) );
        }
    }
}

/**
 * Obtains the values of all variables in the AVMVN proposal on their original scale.
 */
void AVMVNProposal::getAVMVNMemberVariableValuesUntransformed( std::vector<double> *x )
{
    for (size_t i = 0; i < noTransformScalarVariables.size(); ++i)
    {
        x->push_back( noTransformScalarVariables[i]->getValue() );
    }

    for (size_t i = 0; i < logTransformScalarVariables.size(); ++i)
    {
        x->push_back( logTransformScalarVariables[i]->getValue() );
    }

    for (size_t i = 0; i < logitTransformScalarVariables.size(); ++i)
    {
        x->push_back( logitTransformScalarVariables[i]->getValue() );

    }

    for (size_t i = 0; i < noTransformVectorVariables.size(); ++i)
    {

        std::vector<double> tmp = noTransformVectorVariables[i]->getValue();

        for (size_t j=0; j < tmp.size(); ++j)
        {
            x->push_back( tmp[j] );
        }
    }

    for (size_t i = 0; i < logConstrainedSumTransformVectorVariables.size(); ++i)
    {
        std::vector<double> tmp = logConstrainedSumTransformVectorVariables[i]->getValue();

        for (size_t j=0; j < tmp.size(); ++j)
        {
            x->push_back( tmp[j] );
        }
    }
}

/**
 * Resets the values of all variables in the AVMVN proposal.
 */
void AVMVNProposal::resetAVMVNMemberVariableValues( std::vector<double> x )
{
    // Calculate hastings ratio and give variables their new values
    size_t index = 0;
    for (size_t i = 0; i < noTransformScalarVariables.size(); ++i)
    {
        noTransformScalarVariables[i]->getValue() = x[index];
        ++index;
    }

    for (size_t i = 0; i < logTransformScalarVariables.size(); ++i)
    {
        logTransformScalarVariables[i]->getValue() = x[index];
        ++index;
    }

    for (size_t i = 0; i < logitTransformScalarVariables.size(); ++i)
    {
        logitTransformScalarVariables[i]->getValue() = x[index];
        ++index;
    }

    for (size_t i = 0; i < noTransformVectorVariables.size(); ++i)
    {
        std::vector<double> tmp;

        for (size_t j=0; j < noTransformVectorVariables[i]->getValue().size(); ++j)
        {
            tmp.push_back( x[index] );
            ++index;
        }
        noTransformVectorVariables[i]->getValue() = tmp;
    }

    for (size_t i = 0; i < logConstrainedSumTransformVectorVariables.size(); ++i)
    {
        std::vector<double> tmp;

        for (size_t j=0; j < logConstrainedSumTransformVectorVariables[i]->getValue().size(); ++j)
        {
            tmp.push_back( x[index] );
            ++index;
        }
        logConstrainedSumTransformVectorVariables[i]->getValue() = tmp;

    }
}

/**
 * Draw from MVN given that we already have Cholesky decomposition
 *
 * \return The draws.
 */
std::vector<double> AVMVNProposal::rMVNCholesky( std::vector<double> mu, MatrixReal L, RandomNumberGenerator& rng, double scale )
{
    double sqrtScale = sqrt(scale);
    size_t dimension = L.getDim();

    MatrixReal W(dimension, 1, 0.0);
    for (size_t i = 0; i < dimension; ++i)
    {
        W[i][0] = RbStatistics::Normal::rv(0, sqrtScale, rng);
    }

    MatrixReal V = L * W;
    std::vector<double> v = std::vector<double>(dim, 0.0);
    for (size_t i = 0; i < dimension; ++i)
    {
        v[i] = mu[i] + V[i][0];
    }

    return v;

}


/**
 * Perform the proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and Slides the current vale by a scaling factor
 * sf = exp( sigma * u )
 * where sigma is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
double AVMVNProposal::doProposal( void )
{
    ++nTried;

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    lnHastingsratio = 0.0;

    // Get current values in transformed space
    // All variables are accessed in the order they appear in the relevant vector of transformed variables
    // These in turn are accessed in the order: scalar (no transform, log, logit), vector (no transform, sum-constrained-log)
    std::vector<double> x;
    getAVMVNMemberVariableValues(&x);

    std::vector<double> x_untransformed;
    getAVMVNMemberVariableValuesUntransformed(&x_untransformed);

    if ( nTried == 1 )
    // First time using move, setting up components
    {
        dim = x.size();
        MatrixReal vcv(dim);
        for ( size_t i=0; i<dim; ++i )
        {
            storedValues.push_back(x[i]);
            storedValuesUntransformed.push_back(x_untransformed[i]);
            x_bar.push_back(0.0);

            for (size_t j=0; j < dim; ++j)
            {
                vcv[i][j] = 0.0;
            }
            vcv[i][i] = sigma;
        }

        C_emp = vcv;

        vcv.setCholesky(true);
        CholeskyDecomposition& cd = vcv.getCholeskyDecomposition();
        AVMVN_cholesky_L = cd.getLowerCholeskyFactor();
    }
    else
    // Update empirical covariance matrix and averages
    // However, we only change the matrix being used when we tune the variance parameter
    {
        // Store values
        for ( size_t i=0; i<dim; ++i )
        {
            storedValues[i] = x[i];
            storedValuesUntransformed[i] = x_untransformed[i];
        }

        if ( nTried > waitBeforeLearning && updates <= maxUpdates)
        {
            ++updates;

            double n = double(updates);

            for ( size_t i=0; i<dim; ++i )
            {
                // Update covariances first (to save us tracking current and previous averages)
                for (size_t j=i; j < dim; ++j)
                {
                    C_emp[i][j] = 1/n * ( C_emp[i][j] * (n - 1.0) + ((n - 1.0)/n) * ((x[i] - x_bar[i]) * (x[j] - x_bar[j])) );
                    C_emp[j][i] = C_emp[i][j];
                }

                // Update averages
                x_bar[i] = 1/n * x[i] + (n - 1)/n * x_bar[i];
            }

        }

    }

    // Move
    std::vector<double> x_prime = rMVNCholesky(x, AVMVN_cholesky_L, *rng, sigma);

    // This also sets all x to x_prime
    calculateHastingsRatio(x_prime, x);

    return lnHastingsratio;

}


/**
 * Prepare the proposal, e.g., pick the UpDown that we want to change.
 * Here we do not need to do any preparation.
 */
void AVMVNProposal::prepareProposal( void )
{

}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
void AVMVNProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

    o << "sigma = ";
    if (name_only == false)
    {
        o << sigma;
    }

}


/**
 * Remove a variable.
 */
void AVMVNProposal::removeUntransformedScalar( StochasticNode<double> *v )
{

    // add it to the nodes vector
    removeNode( v );

    noTransformScalarVariables.erase(std::remove(noTransformScalarVariables.begin(), noTransformScalarVariables.end(), v), noTransformScalarVariables.end());
}

void AVMVNProposal::removeLogScalar( StochasticNode<double> *v )
{

    // add it to the nodes vector
    removeNode( v );

    logTransformScalarVariables.erase(std::remove(logTransformScalarVariables.begin(), logTransformScalarVariables.end(), v), logTransformScalarVariables.end());
}

void AVMVNProposal::removeLogitScalar( ContinuousStochasticNode *v )
{

    // add it to the nodes vector
    removeNode( v );

    logitTransformScalarVariables.erase(std::remove(logitTransformScalarVariables.begin(), logitTransformScalarVariables.end(), v), logitTransformScalarVariables.end());
}


/**
 * Remove a variable.
 */
void AVMVNProposal::removeUntransformedVector( StochasticNode<RbVector<double> > *v)
{

    // add it to the nodes vector
    removeNode( v );

    noTransformVectorVariables.erase(std::remove(noTransformVectorVariables.begin(), noTransformVectorVariables.end(), v), noTransformVectorVariables.end());
}

void AVMVNProposal::removeLogConstrainedSumVector( StochasticNode<Simplex> *v)
{
    // add it to the nodes vector
    removeNode( v );

    logConstrainedSumTransformVectorVariables.erase(std::remove(logConstrainedSumTransformVectorVariables.begin(), logConstrainedSumTransformVectorVariables.end(), v), logConstrainedSumTransformVectorVariables.end());

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
void AVMVNProposal::undoProposal( void )
{
    resetAVMVNMemberVariableValues(storedValuesUntransformed);

}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
void AVMVNProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    for (size_t i = 0; i < noTransformScalarVariables.size(); ++i)
    {
        if ( noTransformScalarVariables[i] == oldN )
        {
            noTransformScalarVariables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }

    for (size_t i = 0; i < logTransformScalarVariables.size(); ++i)
    {
        if ( logTransformScalarVariables[i] == oldN )
        {
            logTransformScalarVariables[i] = static_cast<StochasticNode<double> *>(newN);
        }
    }

    for (size_t i = 0; i < logitTransformScalarVariables.size(); ++i)
    {
        if ( logitTransformScalarVariables[i] == oldN )
        {
            logitTransformScalarVariables[i] = static_cast<ContinuousStochasticNode *>(newN);
        }
    }

    for (size_t i = 0; i < noTransformVectorVariables.size(); ++i)
    {
        if ( noTransformVectorVariables[i] == oldN )
        {
            noTransformVectorVariables[i] = static_cast<StochasticNode<RbVector<double> > *>(newN);
        }
    }

    for (size_t i = 0; i < logConstrainedSumTransformVectorVariables.size(); ++i)
    {
        if ( logConstrainedSumTransformVectorVariables[i] == oldN )
        {
            logConstrainedSumTransformVectorVariables[i] = static_cast<StochasticNode<Simplex> *>(newN);
        }
    }

}


void AVMVNProposal::setProposalTuningParameter(double tp)
{
    sigma = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
void AVMVNProposal::tune( double rate )
{

    // Update proposal variance
    if ( rate > 0.234 )
    {
        sigma *= (1.0 + ((rate-0.234)/0.764) );
    }
    else
    {
        sigma /= (2.0 - rate/0.234 );
    }

    if ( sigma > 1 ) {
        sigma = fmin(10000, sigma);
    } else {
        sigma = fmax(1/10000, sigma);
    }

    // Update our move's variance-covariance matrix too, then store the Cholesky decomposition
    MatrixReal vcv( dim );

    for (size_t i=0; i<dim; ++i)
    {
        vcv[i][i] = 1.0;
    }

    if ( nTried > waitBeforeUsing )
    {
        for (size_t i=0; i<dim; ++i)
        {
            for (size_t j=i; j<dim; ++j)
            {
                vcv[i][j] = 1.0/dim * sigma * ((1.0 - epsilon) * C_emp[i][j] + epsilon * vcv[i][j]);
                vcv[j][i] = vcv[i][j];
            }
        }
    }
    else
    {
        for (size_t i=0; i<dim; ++i)
        {
            vcv[i][i] = 1.0/dim * sigma;
        }
    }

    vcv.setCholesky(true);
    CholeskyDecomposition& cd = vcv.getCholeskyDecomposition();
    AVMVN_cholesky_L = cd.getLowerCholeskyFactor();

}
