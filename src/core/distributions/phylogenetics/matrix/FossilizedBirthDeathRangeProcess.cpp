#include "FossilizedBirthDeathRangeProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "DistributionExponential.h"
#include "MatrixReal.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathFunctions.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "TimeInterval.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevBayesCore;

/**
 * Constructor. 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    c              Fossil observation counts.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 * \param[in]    re             Augmented age resampling weight.
 */
FossilizedBirthDeathRangeProcess::FossilizedBirthDeathRangeProcess(const DagNode *inspeciation,
                                                                     const DagNode *inextinction,
                                                                     const DagNode *inpsi,
                                                                     const TypedDagNode<double> *inrho,
                                                                     const TypedDagNode< RbVector<double> > *intimes,
                                                                     const std::string &incondition,
                                                                     const std::vector<Taxon> &intaxa,
                                                                     bool complete,
                                                                     bool resample,
                                                                     bool use_bds) :
    TypedDistribution<MatrixReal>(new MatrixReal(intaxa.size(), 2)),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete, resample),
    bds(use_bds)
{
    dirty_gamma = std::vector<bool>(taxa.size(), true);
    gamma_i     = std::vector<size_t>(taxa.size(), 0);
    gamma_links = std::vector<std::vector<bool> >(taxa.size(), std::vector<bool>(taxa.size(), false));

    for(std::vector<const DagNode*>::iterator it = range_parameters.begin(); it != range_parameters.end(); it++)
    {
        addParameter(*it);
    }

    redrawValue();
    updateGamma(true);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
FossilizedBirthDeathRangeProcess* FossilizedBirthDeathRangeProcess::clone( void ) const
{
    return new FossilizedBirthDeathRangeProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathRangeProcess::computeLnProbability( void )
{
    double lnProb = 0.0;

    // prepare the probability computation
    if ( bds == true )
    {
        lnProb = computeLnProbabilityBDS();
    }
    else
    {
        updateGamma();

        lnProb = computeLnProbabilityRanges();
    }

    for( size_t i = 0; i < taxa.size(); i++ )
    {
        // multiply by the number of possible birth locations
        lnProb += log( gamma_i[i] == 0 ? 1 : gamma_i[i] );
    }

    return lnProb;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathRangeProcess::computeLnProbabilityBDS()
{
    // prepare the probability computation
    prepareProbComputation();

    // variable declarations and initialization
    double lnProb = 0.0;

    size_t num_rho_sampled = 0;
    size_t num_rho_unsampled = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        double b = this->getValue()[i][0];
        double d = this->getValue()[i][1];

        double max_age = taxa[i].getMaxAge();
        double min_age = taxa[i].getMinAge();

        double present = times.front();

        // check model constraints
        if ( !( b > o_i[i] && b > d && y_i[i] >= d && d >= present ) )
        {
            return RbConstants::Double::neginf;
        }
        if ( (d > present) != taxa[i].isExtinct() )
        {
            return RbConstants::Double::neginf;
        }

        // count the number of rho-sampled tips
        num_rho_sampled   += (d == present && min_age == present);
        num_rho_unsampled += (d == present && min_age > present);

        if ( dirty_taxa[i] == true )
        {
            size_t bi = findIndex(b);
            size_t di = findIndex(d);

            partial_likelihood[i] = 0.0;

            // include speciation density
            partial_likelihood[i] += log( birth[bi] );

            // skip the rest for extant taxa with no fossil samples
            if ( max_age == present )
            {
                continue;
            }

            // include extinction density
            if (d > present) partial_likelihood[i] += log( death[di] );

            double psi_b_d = 0.0;

            // include poisson density
            for ( size_t j = di; j <= bi; j++ )
            {
                double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

                double dt = std::min(b, t_0) - std::max(d, times[j]);

                partial_likelihood[i] -= (birth[j] + death[j])*dt;

                psi_b_d += fossil[j]*dt;
            }

            // include sampling density
            if ( dirty_psi[i] )
            {
                std::map<TimeInterval, size_t> ages = taxa[i].getOccurrences();

                // if there is a range of fossil ages
                if ( min_age != max_age )
                {
                    double psi_y_o = 0.0;

                    std::vector<double> psi(ages.size(), 0.0);

                    for (size_t j = 0; j < num_intervals; j++)
                    {
                        double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

                        if ( t_0 <= std::max(d,min_age) )
                        {
                            continue;
                        }
                        if ( times[j] >= std::min(b,max_age) )
                        {
                            break;
                        }

                        // increase incomplete sampling psi
                        double dt = std::min(std::min(b,max_age), t_0) - std::max(std::max(d,min_age), times[j]);

                        psi_y_o += fossil[j]*dt;

                        size_t k = 0;
                        // increase running psi total for each observation
                        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                        {
                            if ( Fi->first.getMin() < t_0 && Fi->first.getMax() > times[j] )
                            {
                                double dt = 1.0;

                                // only compute dt if this is a non-singleton
                                if ( Fi->first.getMin() != Fi->first.getMax() )
                                {
                                    dt = std::min(std::min(Fi->first.getMax(), b), t_0) - std::max(std::max(Fi->first.getMin(), d), times[j]);
                                }

                                psi[k] += fossil[j] * dt;
                            }
                        }
                    }

                    // recompute Psi product
                    Psi[i] = 0.0;

                    double count = 0;

                    size_t k = 0;
                    // compute product of psi totals
                    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                    {
                        count += Fi->second;

                        Psi[i] += log(psi[k]) * Fi->second;
                    }

                    if ( complete == true )
                    {
                        // compute poisson density for count
                        Psi[i] -= RbMath::lnFactorial(count);
                        Psi[i] -= psi_b_d;
                    }
                    else
                    {
                        // compute poisson density for count + kappa, kappa >= 0
                        Psi[i] += psi_y_o - psi_b_d;
                        Psi[i] -= count*log(psi_y_o);
                        Psi[i] += log(RbMath::incompleteGamma(psi_y_o, count, true, true));
                    }
                }
                // only one fossil age
                else
                {
                    // include instantaneous sampling density
                    Psi[i] = ages.begin()->second * log(fossil[findIndex(min_age)]);
                }
            }

            if ( condition == "sampling" )
            {
                partial_likelihood[i] -= log(-expm1(-psi_b_d));
            }

            partial_likelihood[i] += Psi[i];
        }

        lnProb += partial_likelihood[i];
    }

    // add the sampled extant tip age term
    if ( homogeneous_rho->getValue() > 0.0)
    {
        lnProb += num_rho_sampled * log( homogeneous_rho->getValue() );
    }
    // add the unsampled extant tip age term
    if ( homogeneous_rho->getValue() < 1.0)
    {
        lnProb += num_rho_unsampled * log( 1.0 - homogeneous_rho->getValue() );
    }

    if ( RbMath::isFinite(lnProb) == false )
    {
        return RbConstants::Double::neginf;
    }

    return lnProb;
}


/**
 * Compute the number of ranges that intersect with range i
 *
 * \param[in]    i      index of range for which to compute gamma
 *
 * \return Small gamma
 */
void FossilizedBirthDeathRangeProcess::updateGamma(bool force)
{
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( dirty_gamma[i] || force )
        {
            double bi = (*this->value)[i][0];
            double di = (*this->value)[i][1];

            if ( force == true ) gamma_i[i] = 0;

            for (size_t j = 0; j < taxa.size(); j++)
            {
                if (i == j) continue;

                double bj = (*this->value)[j][0];
                double dj = (*this->value)[j][1];

                bool linki = ( bi < bj && bi > dj );
                bool linkj = ( bj < bi && bj > di );

                if ( gamma_links[i][j] != linki && force == false )
                {
                    gamma_i[i] += linki ? 1 : -1;
                }
                if ( gamma_links[j][i] != linkj && force == false )
                {
                    gamma_i[j] += linkj ? 1 : -1;
                }

                if ( force == true ) gamma_i[i] += linki;

                gamma_links[i][j] = linki;
                gamma_links[j][i] = linkj;
            }
        }
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
void FossilizedBirthDeathRangeProcess::updateStartEndTimes( void )
{
    origin = 0;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        b_i[i] = (*this->value)[i][0];
        d_i[i] = (*this->value)[i][1];

        origin = std::max(origin, b_i[i]);
    }
}


/**
 * Simulate new speciation times.
 */
void FossilizedBirthDeathRangeProcess::redrawValue(void)
{
    // incorrect placeholder
    // simulation conditioned on the oldest occurrence
    // would require a monte carlo method
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    double max = 0;
    // get the max first occurence
    for (size_t i = 0; i < taxa.size(); i++)
    {
        double o = taxa[i].getMaxAge();
        if ( o > max ) max = o;
    }
    
    max *= 1.1;
    
    if (max == 0.0)
    {
        max = 1.0;
    }

    double present = times.front();

    // get random uniform draws
    for (size_t i = 0; i < taxa.size(); i++)
    {
        // resample oldest occurrence
        resampleAge(i);

        // death time is younger than oldest occurrence and youngest maximum
        double d = taxa[i].isExtinct() ? rng->uniform01()*(std::min(y_i[i], age[i]) - present) + present : present;
        // birth time is older than oldest occurrence
        double b = age[i] + rng->uniform01()*(max - age[i]) + age[i];

        // set values
        (*this->value)[i][0] = b;
        (*this->value)[i][1] = d;
    }
}


void FossilizedBirthDeathRangeProcess::keepSpecialization(const DagNode *toucher)
{
    dirty_gamma = std::vector<bool>(taxa.size(), false);

    AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(toucher);
}

void FossilizedBirthDeathRangeProcess::restoreSpecialization(const DagNode *toucher)
{
    AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathRangeProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    if ( toucher == dag_node )
    {
        if ( touched == false )
        {
            stored_likelihood = partial_likelihood;
            stored_Psi = Psi;

            std::set<size_t> touched_indices = dag_node->getTouchedElementIndices();

            for ( std::set<size_t>::iterator it = touched_indices.begin(); it != touched_indices.end(); it++)
            {
                size_t i = (*it) / taxa.size();

                dirty_gamma[i] = true;
                dirty_psi[i]   = true;
                dirty_taxa[i]  = true;

                if ( resampling == true && resampled == false )
                {
                    resampleAge(i);
                }
            }

        }

        touched = true;
    }
    else
    {
        AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(toucher, touchAll);
    }
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FossilizedBirthDeathRangeProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(oldP, newP);
}
