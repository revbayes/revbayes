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
#include "RbMathLogic.h"
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
                                                                     const std::string &reporting,
                                                                     size_t truncate_at,
                                                                     bool resample,
                                                                     bool use_bds,
                                                                     const TypedDagNode<double> *inorigin,
                                                                     bool report_int) :
    TypedDistribution<MatrixReal>(new MatrixReal(intaxa.size(), 2)),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, reporting, truncate_at, resample, inorigin),
    bds(use_bds)
{
    report_internally = report_int;

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
 * Set the matrix value (e.g. clamping to fixed birth/death ages). A clamp replaces the
 * b/d that redrawValue drew the augmented ages against, so re-clip any age now out of
 * range -- otherwise a clamped chain can start at lnProb = -inf. MCMC moves edit the
 * value in place instead, leaving out-of-range ages for the constraints to reject.
 */
void FossilizedBirthDeathRangeProcess::setValue(MatrixReal *v, bool force)
{
    TypedDistribution<MatrixReal>::setValue(v, force);

    for (size_t i = 0; i < taxa.size(); i++)
    {
        double d  = (*this->value)[i][1];
        double b  = (*this->value)[i][0];

        // oldest age: valid range [max(o_i,d), min(max_age,b))
        double lo = std::max( o_i[i], d );
        double hi = std::min( taxa[i].getMaxAge(), b );
        if ( hi > lo && ( first[i] < lo || first[i] >= b ) )
        {
            first[i] = 0.5 * ( lo + hi );
        }

        // youngest age: valid range [max(d,min_age), min(first,y_i)]
        double lo_y = std::max( d, taxa[i].getMinAge() );
        double hi_y = std::min( first[i], y_i[i] );
        if ( hi_y > lo_y && ( last[i] < lo_y || last[i] > hi_y ) )
        {
            last[i] = 0.5 * ( lo_y + hi_y );
        }
    }
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
        // BDS (Silvestro et al. 2019) treats lineages as independent: no coexistence
        // (gamma) factor. gamma_i is not refreshed here either, so adding one would
        // inject a stale, seed-dependent term.
        lnProb = computeLnProbabilityBDS();
    }
    else
    {
        updateGamma();

        lnProb = computeLnProbabilityRanges();

        for( size_t i = 0; i < taxa.size(); i++ )
        {
            // multiply by the number of possible birth locations
            lnProb += log( gamma_i[i] == 0 ? 1 : gamma_i[i] );
        }
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

        // count the number of rho-sampled tips (see computeLnProbabilityRanges)
        num_rho_sampled += (d == present);

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
                    if ( truncated[i] )
                    {
                    // Truly-exchangeable (uniform subset): the true oldest (tau1 = first[i]) may
                    // be unobserved up to the birth, and interior specimens in (d, tau1) have
                    // their number and ages marginalized. +Lambda because the lineage-wide
                    // e^{-psi_b_d} below already carries the Poisson normalization over [d,b].
                    double o  = first[i];
                    size_t oi = findIndex(o);
                    double Lambda = 0.0;
                    std::vector<double> psi(ages.size(), 0.0);

                    for (size_t j = 0; j < num_intervals; j++)
                    {
                        double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );
                        if ( t_0 <= std::max(d,min_age) ) continue;
                        if ( times[j] >= o ) break;

                        double dL = std::min(o, t_0) - std::max(d, times[j]);
                        if ( dL > 0.0 ) Lambda += fossil[j]*dL;

                        size_t k = 0;
                        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                        {
                            if ( Fi->first.getMin() < t_0 && Fi->first.getMax() > times[j] )
                            {
                                double dt = 1.0;
                                if ( Fi->first.getMin() != Fi->first.getMax() )
                                {
                                    dt = std::min(std::min(Fi->first.getMax(), o), t_0) - std::max(std::max(Fi->first.getMin(), d), times[j]);
                                }
                                psi[k] += fossil[j] * dt;
                            }
                        }
                    }

                    Psi[i] = log(fossil[oi]);
                    int count = 0;
                    double recip_old = 0.0;
                    size_t k = 0;
                    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                    {
                        count += Fi->second;
                        if ( Fi->first.getMin() <= o && Fi->first.getMax() >= o ) recip_old += Fi->second / psi[k];
                        Psi[i] += log(psi[k]) * Fi->second;
                    }

                    double Pk = 0.0;
                    double pmf_count = exp( -Lambda + count*log(Lambda) - RbMath::lnFactorial(count) );
                    double tt = pmf_count;
                    for ( int n = count; n < count + 100000; n++ )
                    {
                        Pk += tt;
                        tt *= Lambda / double(n+1);
                        if ( tt < 1e-17 * Pk && n > (int)Lambda ) break;
                    }
                    double Pk1 = Pk - pmf_count;
                    double bracket = Pk*recip_old + Pk - (double(count)/Lambda)*Pk1;
                    if ( bracket <= 0.0 ) return RbConstants::Double::neginf;
                    Psi[i] += RbMath::lnFactorial(count) - count*log(Lambda) + log(bracket) + Lambda;
                    }
                    else
                    {
                    double psi_int = 0.0;                            // interior sampling rate over (last, first)

                    double o  = first[i];                            // oldest augmented age
                    size_t oi = findIndex(o);
                    double y  = last[i];                             // youngest augmented age
                    size_t yi = findIndex(y);

                    std::vector<double> psi(ages.size(), 0.0);

                    for (size_t j = 0; j < num_intervals; j++)
                    {
                        double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

                        if ( t_0 <= std::max(d,min_age) )
                        {
                            continue;
                        }
                        if ( times[j] >= o )
                        {
                            break;
                        }

                        // the kappa interior spans (last, first) only
                        double dti = std::min(o, t_0) - std::max(y, times[j]);
                        if ( dti > 0.0 ) psi_int += fossil[j]*dti;

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
                                    // observed occurrence sampling rate over its interval (below the oldest)
                                    dt = std::min(std::min(Fi->first.getMax(), o), t_0) - std::max(std::max(Fi->first.getMin(), d), times[j]);
                                }

                                psi[k] += fossil[j] * dt;
                            }
                        }
                    }

                    // include instantaneous sampling density (oldest)
                    Psi[i] = log(fossil[oi]);

                    int count = 0;
                    double recip = 0.0;
                    double recip_young = 0.0;
                    double diag = 0.0;

                    size_t k = 0;
                    // compute factors of the sum over each possible oldest/youngest observation
                    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                    {
                        count += Fi->second;

                        bool eligible_oldest = ( Fi->first.getMax() >= o );
                        bool eligible_youngest = ( Fi->first.getMin() <= y );

                        // compute sum of reciprocal oldest ranges
                        if ( eligible_oldest )
                        {
                            recip += Fi->second / psi[k];
                        }

                        // compute sum of reciprocal youngest ranges
                        if ( eligible_youngest )
                        {
                            recip_young += Fi->second / psi[k];
                        }

                        // intervals that could supply both extremes contribute the diagonal term
                        if ( eligible_oldest && eligible_youngest )
                        {
                            diag += Fi->second / (psi[k]*psi[k]);
                        }

                        // compute product of ranges
                        Psi[i] += log(psi[k]) * Fi->second;
                    }

                    // sum over each possible oldest observation
                    Psi[i] += log(recip);

                    if ( reporting != "firstlast" )
                    {
                        // compute poisson density for count
                        Psi[i] -= RbMath::lnFactorial(count);
                    }
                    else if ( count >= 2 )
                    {
                        // Condition on the oldest and youngest occurrences. Sum over which
                        // observation is the youngest (recip_young) and marginalize the
                        // kappa >= 0 unobserved interior specimens (ages in (last,first)) by
                        // direct summation.
                        if ( !( o >= last[i] && last[i] >= d && last[i] <= y_i[i] && last[i] >= min_age ) )
                        {
                            return RbConstants::Double::neginf;
                        }
                        // youngest instantaneous density + sum over which observation is the youngest,
                        // excluding the diagonal where a single occurrence supplies both extremes
                        Psi[i] += log(fossil[yi]) + log(recip_young - diag/recip);

                        double S1 = 0.0, f = 1.0;
                        for ( size_t kap = 0; kap < 200; kap++ )
                        {
                            S1 += f;
                            f *= psi_int / double(count - 1 + kap);
                            if ( f < 1e-16 * S1 ) break;
                        }
                        Psi[i] += log(S1) - RbMath::lnFactorial(count);
                    }
                    // count == 1 (single occurrence, first == last): no interior term
                    }

                    // no other fossils over the lineage range [d, b]
                    Psi[i] -= psi_b_d;
                }
                // only one fossil age (first == last): include the instantaneous
                // sampling density, the count factorial, and the Poisson
                // normalization over the lineage range.
                else
                {
                    double count = ages.begin()->second;
                    Psi[i]  = count * log(fossil[findIndex(min_age)]);
                    Psi[i] -= RbMath::lnFactorial(count);
                    Psi[i] -= psi_b_d;
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
    double max_birth = 0;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        b_i[i] = (*this->value)[i][0];
        d_i[i] = (*this->value)[i][1];

        max_birth = std::max(max_birth, b_i[i]);
    }

    if ( origin_age != NULL )
    {
        origin = origin_age->getValue();
    }
    else
    {
        origin = max_birth;
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
    // get the max age
    for (size_t i = 0; i < taxa.size(); i++)
    {
        double o = taxa[i].getMaxAge();

        // an unbounded oldest occurrence (max_age = Inf) would put every initial birth time at
        // infinity, so bracket it by the oldest lower bound this taxon does have
        if ( RbMath::isFinite(o) == false )
        {
            o = 0.0;
            const std::map<TimeInterval, size_t>& ages = taxa[i].getOccurrences();
            for ( std::map<TimeInterval, size_t>::const_iterator Fi = ages.begin(); Fi != ages.end(); Fi++ )
            {
                o = std::max( o, Fi->first.getMin() );
            }
        }

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
        // Seed b past every occurrence and d at the present, leaving both clips inert for
        // this call only; updateStartEndTimes overwrites them once the matrix is set.
        b_i[i] = max;
        d_i[i] = present;

        resampleFirstLast(i);

        // the ages are drawn first so d can be placed under them, as the density needs
        // d <= tau_K <= tau_1
        double youngest = std::min( occurrence_counts[i] >= 2 ? last[i] : first[i], y_i[i] );
        double d = taxa[i].isExtinct() ? rng->uniform01()*(youngest - present) + present : present;

        d_i[i] = d;

        // birth time is older than oldest occurrence
        double b = first[i] + rng->uniform01()*(max - first[i]);

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
                size_t i = (*it) / 2; // (birth,death) matrix is N x 2, row-major: linear index / 2 = taxon

                dirty_gamma[i] = true;
                dirty_psi[i]   = true;
                dirty_taxa[i]  = true;

                if ( resampling == true && resampled == false )
                {
                    resampleFirstLast(i);
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
