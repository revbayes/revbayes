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
#include "StochasticNode.h"
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
    condition(incond),
    homogeneous_rho(inrho),
    timeline( intimes ),
    taxa(intaxa),
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

    // cast the pointers from their input parameters
    heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda = dynamic_cast<const TypedDagNode<double >*>(inspeciation);
    heterogeneous_mu = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
    homogeneous_mu = dynamic_cast<const TypedDagNode<double >*>(inextinction);
    heterogeneous_psi = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inpsi);
    homogeneous_psi = dynamic_cast<const TypedDagNode<double >*>(inpsi);

    addParameter( timeline );
    addParameter( homogeneous_rho );
    addParameter( homogeneous_lambda );
    addParameter( heterogeneous_lambda );
    addParameter( homogeneous_mu );
    addParameter( heterogeneous_mu );
    addParameter( homogeneous_psi );
    addParameter( heterogeneous_psi );

    // setup the timeline
    num_intervals = timeline == NULL ? 1 : timeline->getValue().size();

    if ( num_intervals > 1 )
    {
        std::vector<double> times = timeline->getValue();
        std::vector<double> times_sorted_ascending = times;

        sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );

        if ( times != times_sorted_ascending )
        {
            throw(RbException("Interval times must be provided in ascending order"));
        }
    }

    RbException no_timeline_err = RbException("No time intervals provided for heterogeneous fossilized birth death process");

    RbException inconsistent_rates = RbException("Inconsistent number of rates in fossilized birth death process.");

    size_t num_rates = 0;

    if( heterogeneous_lambda != NULL )
    {
        if ( timeline == NULL ) throw(no_timeline_err);

        num_rates = heterogeneous_lambda->getValue().size();
    }
    if( heterogeneous_mu != NULL )
    {
        if ( timeline == NULL ) throw(no_timeline_err);

        if ( num_rates == 0 ) num_rates = heterogeneous_mu->getValue().size();

        if ( heterogeneous_mu->getValue().size() != num_rates ) throw(inconsistent_rates);
    }
    if( heterogeneous_psi != NULL )
    {
        if ( timeline == NULL ) throw(no_timeline_err);

        if ( num_rates == 0 ) num_rates = heterogeneous_psi->getValue().size();

        if ( heterogeneous_psi->getValue().size() != num_rates ) throw(inconsistent_rates);
    }

    if ( num_rates > 0 && num_rates != num_intervals )
    {
        // if all the rate vectors are one longer than the timeline
        // then assume the first time is 0
        if ( num_rates == num_intervals + 1 )
        {
            num_intervals++;
        }
        else
        {
            std::stringstream ss;
            ss << "Number of rates does not match number of time intervals in fossilized birth death process.";
            throw(RbException(ss.str()));
        }
    }

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    o_i = std::vector<double>(taxa.size(), 0.0);
    y_i = std::vector<double>(taxa.size(), RbConstants::Double::inf);

    partial_likelihood = std::vector<double>(taxa.size(), 0.0);

    Psi        = std::vector<double>(taxa.size(), 0.0 );

    dirty_taxa = std::vector<bool>(taxa.size(), true);
    dirty_psi  = std::vector<bool>(taxa.size(), true);

    double max_present = RbConstants::Double::inf;

    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        std::map<TimeInterval, size_t> ages = taxa[i].getOccurrences();
        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++ )
        {
            // find the oldest minimum age
            o_i[i] = std::max(Fi->first.getMin(), o_i[i]);
            // find the youngest maximum age
            y_i[i] = std::min(Fi->first.getMax(), y_i[i]);

            max_present = std::min(max_present, y_i[i]);
        }
    }

    prepareProbComputation();
    
    if ( times.front() > max_present )
    {
        throw(RbException("Timeline start time is older than youngest fossil age."));
    }

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
 * return the index i so that t_i <= t < t_{i+1}
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 = 0.0
 * t_l is origin
 */
size_t BirthDeathRateShiftsProcess::findIndex(double t) const
{
    return std::prev(std::upper_bound( times.begin(), times.end(), t)) - times.begin();
}


/**
 * Simulate new speciation times.
 */
void BirthDeathRateShiftsProcess::redrawValue(void)
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
        double b = taxa[i].getMaxAge() + rng->uniform01()*(max - taxa[i].getMaxAge()) + taxa[i].getMaxAge();
        double d = taxa[i].isExtinct() ? rng->uniform01()*(taxa[i].getMinAge() - present) + present : present;

        (*this->value)[i][0] = b;
        (*this->value)[i][1] = d;
    }
}


/**
 *
 *
 */
void BirthDeathRateShiftsProcess::prepareProbComputation()
{
    if ( homogeneous_lambda != NULL )
    {
        birth = std::vector<double>(num_intervals, homogeneous_lambda->getValue() );
    }
    else
    {
        birth = heterogeneous_lambda->getValue();
    }
    if ( homogeneous_mu != NULL )
    {
        death = std::vector<double>(num_intervals, homogeneous_mu->getValue() );
    }
    else
    {
        death = heterogeneous_mu->getValue();
    }
    if ( homogeneous_psi != NULL )
    {
        fossil = std::vector<double>(num_intervals, homogeneous_psi->getValue() );
    }
    else
    {
        fossil = heterogeneous_psi->getValue();
    }

    if ( timeline != NULL )
    {
        times = timeline->getValue();
    }
    else
    {
        times.clear();
    }

    if ( times.size() < num_intervals )
    {
        times.insert(times.begin(), 0.0);
    }
}


void BirthDeathRateShiftsProcess::keepSpecialization(DagNode *toucher)
{
    dirty_psi  = std::vector<bool>(taxa.size(), false);
    dirty_taxa = std::vector<bool>(taxa.size(), false);

    touched = false;
}


void BirthDeathRateShiftsProcess::restoreSpecialization(DagNode *toucher)
{
    partial_likelihood = stored_likelihood;
    Psi = stored_Psi;

    dirty_psi  = std::vector<bool>(taxa.size(), false);
    dirty_taxa = std::vector<bool>(taxa.size(), false);

    touched = false;
}


void BirthDeathRateShiftsProcess::touchSpecialization(DagNode *toucher, bool touchAll)
{
    if ( touched == false )
    {
        stored_likelihood = partial_likelihood;
        stored_Psi = Psi;

        if ( toucher == timeline || toucher == homogeneous_psi || toucher == heterogeneous_psi || touchAll )
        {
            dirty_taxa = std::vector<bool>(taxa.size(), true);
            dirty_psi  = std::vector<bool>(taxa.size(), true);
        }
        else if ( toucher == dag_node )
        {
            std::set<size_t> touched_indices = dag_node->getTouchedElementIndices();

            for ( std::set<size_t>::iterator it = touched_indices.begin(); it != touched_indices.end(); it++)
            {
                size_t i = (*it) / taxa.size();

                dirty_psi[i]  = true;
                dirty_taxa[i]  = true;
            }
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
