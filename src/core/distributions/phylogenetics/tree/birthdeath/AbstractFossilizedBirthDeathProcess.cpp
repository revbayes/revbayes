#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

#include "AbstractFossilizedBirthDeathProcess.h"
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
AbstractFossilizedBirthDeathProcess::AbstractFossilizedBirthDeathProcess(const DagNode *inspeciation,
                                                                         const DagNode *inextinction,
                                                                         const DagNode *inpsi,
                                                                         const TypedDagNode<double> *inrho,
                                                                         const TypedDagNode< RbVector<double> > *intimes,
                                                                         const std::vector<Taxon> &intaxa,
                                                                         bool c) :
    ascending(true), homogeneous_rho(inrho), timeline( intimes ), fbd_taxa(intaxa), complete(c), origin(0.0), touched(false)
{
    // initialize all the pointers to NULL
    homogeneous_lambda             = NULL;
    homogeneous_mu                 = NULL;
    homogeneous_psi                = NULL;
    heterogeneous_lambda           = NULL;
    heterogeneous_mu               = NULL;
    heterogeneous_psi              = NULL;

    RbException no_timeline_err = RbException("No time intervals provided for piecewise constant fossilized birth death process");

    range_parameters.push_back( homogeneous_rho );

    heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

    range_parameters.push_back( homogeneous_lambda );
    range_parameters.push_back( heterogeneous_lambda );

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

    range_parameters.push_back( homogeneous_mu );
    range_parameters.push_back( heterogeneous_mu );

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

    range_parameters.push_back( homogeneous_psi );
    range_parameters.push_back( heterogeneous_psi );

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

    range_parameters.push_back( timeline );

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

    b_i = std::vector<double>(fbd_taxa.size(), 0.0);
    d_i = std::vector<double>(fbd_taxa.size(), 0.0);

    p_i         = std::vector<double>(num_intervals, 1.0);
    q_i         = std::vector<double>(num_intervals, 0.0);
    q_tilde_i   = std::vector<double>(num_intervals, 0.0);

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    updateIntervals();

    partial_likelihood = std::vector<double>(fbd_taxa.size(), 0.0);

    o_i         = std::vector<double>(fbd_taxa.size(), 0.0);
    y_i         = std::vector<size_t>(fbd_taxa.size(), 0.0);
    x_i         = std::vector<std::vector<double> >(fbd_taxa.size(), std::vector<double>() );
    nu_j        = std::vector<std::vector<double> >(fbd_taxa.size(), std::vector<double>() );
    Psi_i       = std::vector<std::vector<double> >(fbd_taxa.size(), std::vector<double>() );

    dirty_taxa = std::vector<bool>(fbd_taxa.size(), true);
    dirty_psi = std::vector<bool>(fbd_taxa.size(), true);

    for ( size_t i = 0; i < fbd_taxa.size(); i++ )
    {
        std::map<TimeInterval, size_t> ages = fbd_taxa[i].getAges();

        std::set<double> uv;
        double oldest_y = 0.0;

        // get sorted unique uncertainty range ages
        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++ )
        {
            uv.insert(Fi->first.getMin());
            uv.insert(Fi->first.getMax());

            oldest_y = std::max(Fi->first.getMin(), oldest_y);
        }

        // put the sorted ages in a vector
        for ( std::set<double>::iterator it = uv.begin(); it != uv.end(); it++ )
        {
            x_i[i].push_back(*it);
        }

        // compute nu
        for ( size_t j = 0; j < x_i[i].size(); j++ )
        {
            size_t nu = 0;

            for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++ )
            {
                // count the number of ranges with maximum > xj
                nu += ( Fi->first.getMax() > x_i[i][j] ) * Fi->second;
            }

            if ( x_i[i][j] == oldest_y )
            {
                y_i[i] = j;
            }

            nu_j[i].push_back(nu);
        }

        std::vector<double> x;
        std::set_union( times.begin(), times.end(), x_i[i].begin(), x_i[i].end(), std::back_inserter(x) );

        Psi_i[i] = std::vector<double>(x.size(), 0.0);
    }
}

/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double AbstractFossilizedBirthDeathProcess::computeLnProbabilityRanges( bool force )
{
    // prepare the probability computation
    updateIntervals();

    // variable declarations and initialization
    double lnProbTimes = 0.0;

    size_t num_extant_sampled = 0;
    size_t num_extant_unsampled = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < fbd_taxa.size(); ++i)
    {
        double b = b_i[i];
        double d = d_i[i];

        double o = fbd_taxa[i].getMaxAge();
        double y = fbd_taxa[i].getMinAge();

        // check model constraints
        if ( !( b > o && y >= d && d >= 0.0 ) )
        {
            return RbConstants::Double::neginf;
        }
        if ( (d > 0.0) != fbd_taxa[i].isExtinct() )
        {
            return RbConstants::Double::neginf;
        }

        // count the number of rho-sampled tips
        num_extant_sampled   += (d == 0.0 && y == 0.0);  // l
        num_extant_unsampled += (d == 0.0 && y > 0.0); // n - m - l

        if ( dirty_taxa[i] == true || force )
        {
            size_t bi = l(b);
            size_t oi = l(o_i[i]);
            size_t di = l(d);

            partial_likelihood[i] = 0.0;

            // include speciation density
            partial_likelihood[i] += log( birth[bi] );

            // multiply by q at the birth time
            partial_likelihood[i] += q(bi, b);

            // include intermediate q terms
            for (size_t j = bi; j < oi; j++)
            {
                partial_likelihood[i] += q_i[j];
            }

            // replace q terms at oldest occurrence
            partial_likelihood[i] += q(oi, o_i[i], true) - q(oi, o_i[i]);

            // include intermediate q_tilde terms
            for (size_t j = oi; j < di; j++)
            {
                partial_likelihood[i] += q_tilde_i[j];
            }

            // divide by q_tilde at the death time
            partial_likelihood[i] -= q( di, d, true);

            // include extinction density
            if (d > 0.0) partial_likelihood[i] += log( death[di] );


            if ( dirty_psi[i] || force )
            {
                std::map<TimeInterval, size_t> ages = fbd_taxa[i].getAges();

                double psi_y_o = 0.0;

                std::vector<double> psi(ages.size(), 0.0);

                for (size_t interval = num_intervals; interval > 0; interval--)
                {
                    size_t j = interval - 1;

                    double t_0 = ( j > 0 ? times[j-1] : RbConstants::Double::inf );

                    if ( t_0 <= fbd_taxa[i].getMinAge() )
                    {
                        continue;
                    }
                    if ( times[j] >= o_i[i] )
                    {
                        break;
                    }

                    // increase incomplete sampling psi
                    double dt = std::min(o_i[i], t_0) - std::max(fbd_taxa[i].getMinAge(), times[j]);

                    psi_y_o += fossil[j]*dt;

                    size_t k = 0;
                    // increase running psi total for each observation
                    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                    {
                        if ( Fi->first.getMin() >= t_0 || Fi->first.getMax() <= times[j] )
                        {
                            continue;
                        }

                        dt = std::min(std::min(Fi->first.getMax(), o_i[i]), t_0) - std::max(Fi->first.getMin(), times[j]);

                        psi[k] += fossil[j]*dt;
                    }
                }

                // include sampling density
                Psi_i[i][0] = log(fossil[oi]);

                double recip = 0.0;

                size_t k = 0;
                // compute factors of the sum over each possible oldest observation
                for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                {
                    // compute sum of reciprocal ranges
                    if ( Fi->first.getMax() >= o_i[i] )
                    {
                        recip += Fi->second / psi[k];
                    }

                    // compute product of ranges
                    Psi_i[i][0] += log(psi[k]) * Fi->second;
                }

                // sum over each possible oldest observation
                Psi_i[i][0] += log(recip);

                // multiply by e^psi_y_o
                Psi_i[i][0] += ( complete ? 0.0 : psi_y_o );
            }

            partial_likelihood[i] += Psi_i[i][0];
        }

        lnProbTimes += partial_likelihood[i];
    }

    // (the origin is not a speciation event)
    lnProbTimes -= log( birth[l(origin)] );

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


double AbstractFossilizedBirthDeathProcess::getExtinctionRate( size_t index ) const
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


double AbstractFossilizedBirthDeathProcess::getFossilSamplingRate( size_t index ) const
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


double AbstractFossilizedBirthDeathProcess::getIntervalTime( size_t index ) const
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


double AbstractFossilizedBirthDeathProcess::getSpeciationRate( size_t index ) const
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
size_t AbstractFossilizedBirthDeathProcess::l(double t) const
{
    return times.rend() - std::upper_bound( times.rbegin(), times.rend(), t);
}


/**
 * p_i(t)
 */
double AbstractFossilizedBirthDeathProcess::p( size_t i, double t ) const
{
    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == num_intervals - 1 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];
    
    double diff = b - d - f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*b*f);
    double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

    double ln_e = -A*dt;

    double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);
    
    return (b + d + f - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);
}


/**
 * q_i(t)
 */
double AbstractFossilizedBirthDeathProcess::q( size_t i, double t, bool tilde ) const
{
    
    if ( t == 0.0 ) return 1.0;
    
    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == num_intervals - 1 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];
    
    double diff = b - d - f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*b*f);
    double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

    double ln_e = -A*dt;

    double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

    double q = log(4.0) + ln_e - 2.0*log(tmp);

    if (tilde) q = 0.5 * (q - (b+d+f)*dt);
    
    return q;
}


/**
 *
 *
 */
void AbstractFossilizedBirthDeathProcess::redrawOldestOccurrence(size_t i)
{
    dirty_taxa[i] = true;
    dirty_psi[i]  = true;
    o_i[i] = GLOBAL_RNG->uniform01()*(fbd_taxa[i].getMaxAge() - x_i[i][y_i[i]]) + x_i[i][y_i[i]];
}


void AbstractFossilizedBirthDeathProcess::keepSpecialization(DagNode *toucher)
{
    dirty_psi  = std::vector<bool>(fbd_taxa.size(), false);
    dirty_taxa = std::vector<bool>(fbd_taxa.size(), false);

    touched = false;
}


void AbstractFossilizedBirthDeathProcess::restoreSpecialization(DagNode *toucher)
{
    partial_likelihood = stored_likelihood;
    o_i = stored_o_i;
    Psi_i = stored_Psi_i;

    dirty_psi  = std::vector<bool>(fbd_taxa.size(), false);
    dirty_taxa = std::vector<bool>(fbd_taxa.size(), false);

    touched = false;
}


void AbstractFossilizedBirthDeathProcess::touchSpecialization(DagNode *toucher, bool touchAll)
{
    if ( touched == false )
    {
        stored_likelihood = partial_likelihood;
        stored_o_i = o_i;
        stored_Psi_i = Psi_i;

        dirty_taxa = std::vector<bool>(fbd_taxa.size(), true);

        if ( toucher == timeline || toucher == homogeneous_psi || toucher == heterogeneous_psi || touchAll )
        {
            dirty_psi  = std::vector<bool>(fbd_taxa.size(), true);
        }
    }

    touched = true;
}


/**
 *
 *
 */
void AbstractFossilizedBirthDeathProcess::updateIntervals()
{
    for (size_t interval = num_intervals; interval > 0; interval--)
    {
        size_t i = interval - 1;

        double b = getSpeciationRate(i);
        double d = getExtinctionRate(i);
        double f = getFossilSamplingRate(i);
        double ti = getIntervalTime(i);

        birth[i] = b;
        death[i] = d;
        fossil[i] = f;
        times[i] = ti;

        if (i > 0)
        {
            double r = (i == num_intervals - 1 ? homogeneous_rho->getValue() : 0.0);
            double t = getIntervalTime(i-1);

            double diff = b - d - f;
            double dt   = t - ti;

            double A = sqrt( diff*diff + 4.0*b*f);
            double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

            double ln_e = -A*dt;

            double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

            q_i[i-1]       = log(4.0) + ln_e - 2.0*log(tmp);
            q_tilde_i[i-1] = 0.5 * ( q_i[i-1] - (b+d+f)*dt );
            p_i[i-1]       = (b + d + f - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);
        }
    }
}


/**
 * Swap the parameters held by this distribution.
 * 
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void AbstractFossilizedBirthDeathProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
