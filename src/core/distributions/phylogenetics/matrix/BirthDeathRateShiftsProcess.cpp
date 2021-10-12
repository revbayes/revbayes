#include "BirthDeathRateShiftsProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

#include "DagNode.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RlUserInterface.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TypedDagNode.h"

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
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Taxa.
 */
BirthDeathRateShiftsProcess::BirthDeathRateShiftsProcess(const DagNode *inspeciation,
                                                         const DagNode *inextinction,
                                                         const DagNode *inpsi,
                                                         const TypedDagNode<double> *inrho,
                                                         const TypedDagNode< RbVector<double> > *intimes,
                                                         const std::string& incond,
                                                         const std::vector<Taxon> &intaxa,
                                                         bool c) :
    TypedDistribution<MatrixReal>(new MatrixReal(intaxa.size(), 2)),
    ascending(true),
    condition(incond),
    homogeneous_rho(inrho),
    timeline( intimes ),
    bd_taxa(intaxa),
    complete(c),
    touched(false)
{
    // initialize all the pointers to NULL
    homogeneous_lambda             = NULL;
    homogeneous_mu                 = NULL;
    homogeneous_psi                = NULL;
    heterogeneous_lambda           = NULL;
    heterogeneous_mu               = NULL;
    heterogeneous_psi              = NULL;

    RbException no_timeline_err = RbException("No time intervals provided for time-heterogeneous birth death with rate shifts process");

    addParameter( homogeneous_rho );

    heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

    addParameter( homogeneous_lambda );
    addParameter( heterogeneous_lambda );

    if( heterogeneous_lambda != NULL )
    {
        if( timeline == NULL ) throw(no_timeline_err);

        if (heterogeneous_lambda->getValue().size() != timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of speciation rates (" << heterogeneous_lambda->getValue().size() << ") does not match number of time intervals (" << timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
    }


    heterogeneous_mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
    homogeneous_mu = dynamic_cast<const TypedDagNode<double >*>(inextinction);

    addParameter( homogeneous_mu );
    addParameter( heterogeneous_mu );

    if( heterogeneous_mu != NULL )
    {
        if( timeline == NULL ) throw(no_timeline_err);

        if (heterogeneous_mu->getValue().size() != timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of extinction rates (" << heterogeneous_mu->getValue().size() << ") does not match number of time intervals (" << timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
    }


    heterogeneous_psi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inpsi);
    homogeneous_psi = dynamic_cast<const TypedDagNode<double >*>(inpsi);

    addParameter( homogeneous_psi );
    addParameter( heterogeneous_psi );

    if( heterogeneous_psi != NULL )
    {
        if( timeline == NULL ) throw(no_timeline_err);

        if (heterogeneous_psi->getValue().size() != timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of fossil sampling rates (" << heterogeneous_psi->getValue().size() << ") does not match number of time intervals (" << timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
    }

    addParameter( timeline );

    num_intervals = timeline == NULL ? 1 : timeline->getValue().size()+1;

    if ( num_intervals > 1 )
    {
        std::vector<double> times = timeline->getValue();
        std::vector<double> times_sorted_ascending = times;
        std::vector<double> times_sorted_descending = times;

        sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
        sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

        if( times == times_sorted_descending )
        {
            ascending = false;
        }
        else if ( times != times_sorted_ascending )
        {
            throw(RbException("Interval times must be provided in order"));
        }
    }

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    partial_likelihood = std::vector<double>(bd_taxa.size(), 0.0);

    Psi_i       = std::vector<double>(bd_taxa.size(), 0.0 );

    dirty_taxa = std::vector<bool>(bd_taxa.size(), true);
    dirty_psi = std::vector<bool>(bd_taxa.size(), true);

    updateIntervals();
    redrawValue();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself
 */
BirthDeathRateShiftsProcess* BirthDeathRateShiftsProcess::clone( void ) const
{
    return new BirthDeathRateShiftsProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double BirthDeathRateShiftsProcess::computeLnProbability()
{
    // prepare the probability computation
    updateIntervals();

    // variable declarations and initialization
    double lnProbTimes = 0.0;

    size_t num_extant_sampled = 0;
    size_t num_extant_unsampled = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < bd_taxa.size(); ++i)
    {
        double b = this->getValue()[i][0];
        double d = this->getValue()[i][1];

        double o = bd_taxa[i].getMaxAge();
        double y = bd_taxa[i].getMinAge();

        // check model constraints
        if ( !( b > o && y >= d && d >= 0.0 ) )
        {
            return RbConstants::Double::neginf;
        }
        if ( (d > 0.0) != bd_taxa[i].isExtinct() )
        {
            return RbConstants::Double::neginf;
        }

        // count the number of rho-sampled tips
        num_extant_sampled   += (d == 0.0 && y == 0.0);
        num_extant_unsampled += (d == 0.0 && y > 0.0);

        if ( dirty_taxa[i] == true )
        {
            size_t bi = l(b);
            size_t di = l(d);

            partial_likelihood[i] = 0.0;

            // include speciation density
            partial_likelihood[i] += log( birth[bi] );

            // include extinction density
            if (d > 0.0) partial_likelihood[i] += log( death[di] );

            double psi_b_d = 0.0;

            // include poisson density
            for ( size_t j = bi; j <= di; j++ )
            {
                double t_0 = ( j > 0 ? times[j-1] : RbConstants::Double::inf );

                // increase incomplete sampling psi
                double dt = std::min(b, t_0) - std::max(d, times[j]);

                partial_likelihood[i] -= (birth[j] + death[j])*dt;

                psi_b_d += fossil[j]*dt;
            }

            // include sampling density
            if ( dirty_psi[i] )
            {
                std::map<TimeInterval, size_t> ages = bd_taxa[i].getAges();

                double psi_y_o = 0.0;

                std::vector<double> psi(ages.size(), 0.0);

                for (size_t interval = num_intervals; interval > 0; interval--)
                {
                    size_t j = interval - 1;

                    double t_0 = ( j > 0 ? times[j-1] : RbConstants::Double::inf );

                    if ( t_0 <= bd_taxa[i].getMinAge() )
                    {
                        continue;
                    }
                    if ( times[j] >= bd_taxa[i].getMaxAge() )
                    {
                        break;
                    }

                    // increase incomplete sampling psi
                    double dt = std::min(bd_taxa[i].getMaxAge(), t_0) - std::max(bd_taxa[i].getMinAge(), times[j]);

                    psi_y_o += fossil[j]*dt;

                    size_t k = 0;
                    // increase running psi total for each observation
                    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                    {
                        if ( Fi->first.getMin() >= t_0 || Fi->first.getMax() <= times[j] )
                        {
                            continue;
                        }

                        dt = std::min(Fi->first.getMax(), t_0) - std::max(Fi->first.getMin(), times[j]);

                        psi[k] += fossil[j]*dt;
                    }
                }

                // recompute Psi product
                Psi_i[i] = 0.0;

                double count = 0;

                size_t k = 0;
                // compute product of psi totals
                for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                {
                    count += Fi->second;

                    Psi_i[i] += log(psi[k]) * Fi->second;
                }

                if ( complete == true )
                {
                    // compute poisson density for count
                    Psi_i[i] -= RbMath::lnFactorial(count);
                    Psi_i[i] -= psi_b_d;
                }
                else
                {
                    // compute poisson density for count + kappa, kappa >= 0
                    Psi_i[i] += psi_y_o - psi_b_d;
                    Psi_i[i] -= count*log(psi_b_d);
                    Psi_i[i] += log(RbMath::incompleteGamma(psi_b_d, count, true, true));
                }
            }

            if ( condition == "sampling" )
            {
                partial_likelihood[i] -= log(-expm1(-psi_b_d));
            }

            partial_likelihood[i] += Psi_i[i];
        }

        lnProbTimes += partial_likelihood[i];
    }

    // add the sampled extant tip age term
    if ( homogeneous_rho->getValue() > 0.0)
    {
        lnProbTimes += num_extant_sampled * log( homogeneous_rho->getValue() );
    }
    // add the unsampled extant tip age term
    if ( homogeneous_rho->getValue() < 1.0)
    {
        lnProbTimes += num_extant_unsampled * log( 1.0 - homogeneous_rho->getValue() );
    }

    if ( RbMath::isFinite(lnProbTimes) == false )
    {
        return RbConstants::Double::neginf;
    }

    return lnProbTimes;
}


double BirthDeathRateShiftsProcess::getExtinctionRate( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_mu != NULL )
    {
        return homogeneous_mu->getValue();
    }
    else
    {
        size_t num = heterogeneous_mu->getValue().size();

        if (index >= num)
        {
            throw(RbException("Extinction rate index out of bounds"));
        }
        return ascending ? heterogeneous_mu->getValue()[num - 1 - index] : heterogeneous_mu->getValue()[index];
    }
}


double BirthDeathRateShiftsProcess::getFossilSamplingRate( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_psi != NULL )
    {
        return homogeneous_psi->getValue();
    }
    else
    {
        size_t num = heterogeneous_psi->getValue().size();

        if (index >= num)
        {
            throw(RbException("Fossil sampling rate index out of bounds"));
        }
        return ascending ? heterogeneous_psi->getValue()[num - 1 - index] : heterogeneous_psi->getValue()[index];
    }
}


double BirthDeathRateShiftsProcess::getIntervalTime( size_t index ) const
{

    if ( index == num_intervals - 1 )
    {
        return 0.0;
    }
    // remove the old parameter first
    else if ( timeline != NULL )
    {
        size_t num = timeline->getValue().size();

        if (index >= num)
        {
            throw(RbException("Interval time index out of bounds"));
        }
        return ascending ? timeline->getValue()[num - 1 - index] : timeline->getValue()[index];
    }
    else
    {
        throw(RbException("Interval time index out of bounds"));
    }
}


double BirthDeathRateShiftsProcess::getSpeciationRate( size_t index ) const
{

    // remove the old parameter first
    if ( homogeneous_lambda != NULL )
    {
        return homogeneous_lambda->getValue();
    }
    else
    {
        size_t num = heterogeneous_lambda->getValue().size();

        if (index >= num)
        {
            throw(RbException("Speciation rate index out of bounds"));
        }
        return ascending ? heterogeneous_lambda->getValue()[num - 1 - index] : heterogeneous_lambda->getValue()[index];
    }
}


/**
 * return the index i so that t_{i-1} > t >= t_i
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 is origin
 * t_l = 0.0
 */
size_t BirthDeathRateShiftsProcess::l(double t) const
{
    return times.rend() - std::upper_bound( times.rbegin(), times.rend(), t);
}


/**
 * Simulate new speciation times.
 */
void BirthDeathRateShiftsProcess::redrawValue(void)
{
    // incorrect placeholder

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    double max = 0;
    // get the max first occurence
    for (size_t i = 0; i < bd_taxa.size(); i++)
    {
        double o = bd_taxa[i].getMaxAge();
        if ( o > max ) max = o;
    }

    max *= 1.1;

    if (max == 0.0)
    {
        max = 1.0;
    }

    // get random uniform draws
    for (size_t i = 0; i < bd_taxa.size(); i++)
    {
        double b = bd_taxa[i].getMaxAge() + rng->uniform01()*(max - bd_taxa[i].getMaxAge());
        double d = bd_taxa[i].isExtinct() ? rng->uniform01()*bd_taxa[i].getMinAge() : 0.0;

        (*this->value)[i][0] = b;
        (*this->value)[i][1] = d;
    }
}


/**
 *
 *
 */
void BirthDeathRateShiftsProcess::updateIntervals()
{
    for (size_t interval = num_intervals; interval > 0; interval--)
    {
        size_t i = interval - 1;

        birth[i] = getSpeciationRate(i);
        death[i] = getExtinctionRate(i);
        fossil[i] = getFossilSamplingRate(i);
        times[i] = getIntervalTime(i);
    }
}


void BirthDeathRateShiftsProcess::keepSpecialization(DagNode *toucher)
{
    dirty_psi  = std::vector<bool>(bd_taxa.size(), false);
    dirty_taxa = std::vector<bool>(bd_taxa.size(), false);

    touched = false;
}


void BirthDeathRateShiftsProcess::restoreSpecialization(DagNode *toucher)
{
    partial_likelihood = stored_likelihood;
    Psi_i = stored_Psi_i;

    dirty_psi  = std::vector<bool>(bd_taxa.size(), false);
    dirty_taxa = std::vector<bool>(bd_taxa.size(), false);

    touched = false;
}


void BirthDeathRateShiftsProcess::touchSpecialization(DagNode *toucher, bool touchAll)
{
    if ( touched == false )
    {
        stored_likelihood = partial_likelihood;
        stored_Psi_i = Psi_i;

        dirty_taxa = std::vector<bool>(bd_taxa.size(), true);

        if ( toucher == timeline || toucher == homogeneous_psi || toucher == heterogeneous_psi || touchAll )
        {
            dirty_psi  = std::vector<bool>(bd_taxa.size(), true);
        }
    }

    touched = true;
}


/**
 * Swap the parameters held by this distribution.
 * 
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void BirthDeathRateShiftsProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == heterogeneous_lambda)
    {
        heterogeneous_lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_mu)
    {
        heterogeneous_mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_psi)
    {
        heterogeneous_psi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda)
    {
        homogeneous_lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_mu)
    {
        homogeneous_mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_psi)
    {
        homogeneous_psi = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_rho)
    {
        homogeneous_rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == timeline)
    {
        timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
}
