#include "AbstractFossilizedBirthDeathProcess.h"

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
AbstractFossilizedBirthDeathProcess::AbstractFossilizedBirthDeathProcess(const DagNode *inspeciation,
                                                                         const DagNode *inextinction,
                                                                         const DagNode *inpsi,
                                                                         const TypedDagNode<double> *inrho,
                                                                         const TypedDagNode< RbVector<double> > *intimes,
                                                                         const std::vector<Taxon> &intaxa,
                                                                         bool c,
                                                                         bool ex,
                                                                         double re) :
    fbd_taxa(intaxa),
    homogeneous_rho(inrho),
    timeline( intimes ),
    origin(0.0),
    complete(c),
    extended(ex),
    touched(false),
    resampling(re)
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

    // add the parameters to the model
    range_parameters.push_back( timeline );
    range_parameters.push_back( homogeneous_rho );
    range_parameters.push_back( homogeneous_lambda );
    range_parameters.push_back( heterogeneous_lambda );
    range_parameters.push_back( homogeneous_mu );
    range_parameters.push_back( heterogeneous_mu );
    range_parameters.push_back( homogeneous_psi );
    range_parameters.push_back( heterogeneous_psi );

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

    if ( num_rates != num_intervals )
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

    b_i = std::vector<double>(fbd_taxa.size(), 0.0);
    d_i = std::vector<double>(fbd_taxa.size(), 0.0);
    o_i = std::vector<double>(fbd_taxa.size(), 0.0);
    y_i = std::vector<double>(fbd_taxa.size(), RbConstants::Double::inf);

    p_i         = std::vector<double>(num_intervals, 1.0);
    q_i         = std::vector<double>(num_intervals, 0.0);
    q_tilde_i   = std::vector<double>(num_intervals, 0.0);

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    partial_likelihood = std::vector<double>(fbd_taxa.size(), 0.0);

    age        = std::vector<double>(fbd_taxa.size(), 0.0);
    Psi         = std::vector<double>(fbd_taxa.size(), 0.0 );

    dirty_taxa = std::vector<bool>(fbd_taxa.size(), true);
    dirty_psi = std::vector<bool>(fbd_taxa.size(), true);

    for ( size_t i = 0; i < fbd_taxa.size(); i++ )
    {
        std::map<TimeInterval, size_t> ages = fbd_taxa[i].getAges();
        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++ )
        {
            // find the oldest minimum age
            o_i[i] = std::max(Fi->first.getMin(), o_i[i]);
            // find the youngest maximum age
            y_i[i] = std::min(Fi->first.getMax(), y_i[i]);
        }
    }

    prepareProbComputation();
}

/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double AbstractFossilizedBirthDeathProcess::computeLnProbabilityRanges( bool force )
{
    // prepare the probability computation
    prepareProbComputation();

    // variable declarations and initialization
    double lnProbTimes = 0.0;

    size_t num_extant_sampled = 0;
    size_t num_extant_unsampled = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < fbd_taxa.size(); ++i)
    {
        double b = b_i[i];
        double d = d_i[i];

        double o = age[i];

        double max_age = fbd_taxa[i].getMaxAge();
        double min_age = fbd_taxa[i].getMinAge();

        // check model constraints
        if ( !( b > o && o > o_i[i] && y_i[i] >= d && (( extended && d >= 0.0) || (!extended && d > min_age ) )) )
        {
            return RbConstants::Double::neginf;
        }
        if ( (d > 0.0) != fbd_taxa[i].isExtinct() )
        {
            return RbConstants::Double::neginf;
        }

        // count the number of rho-sampled tips
        num_extant_sampled   += (d == 0.0 && min_age == 0.0);  // l
        num_extant_unsampled += (d == 0.0 && min_age > 0.0); // n - m - l

        if ( dirty_taxa[i] == true || force )
        {
            size_t bi = findIndex(b);
            size_t oi = findIndex(o);
            size_t di = findIndex(d);

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

            // skip the rest for extant taxa with no fossil samples
            if ( max_age == 0.0 )
            {
                lnProbTimes += partial_likelihood[i];

                continue;
            }

            // replace q terms at oldest occurrence
            partial_likelihood[i] += q(oi, o, true) - q(oi, o);

            // include intermediate q_tilde terms
            for (size_t j = oi; j < di; j++)
            {
                partial_likelihood[i] += q_tilde_i[j];
            }

            // divide by q_tilde at the death time
            partial_likelihood[i] -= q( di, d, true);

            // include extinction density
            if ( d > 0.0 ) partial_likelihood[i] += extended ? log( death[di] ) : log( p(di, d) );

            if ( dirty_psi[i] || force )
            {
                std::map<TimeInterval, size_t> ages = fbd_taxa[i].getAges();

                // if there is a range of fossil ages
                if ( min_age != max_age )
                {
                    double psi_y_o = 0.0;

                    std::vector<double> psi(ages.size(), 0.0);

                    for (size_t interval = num_intervals; interval > 0; interval--)
                    {
                        size_t j = interval - 1;

                        double t_0 = ( j > 0 ? times[j-1] : RbConstants::Double::inf );

                        if ( t_0 <= min_age )
                        {
                            continue;
                        }
                        if ( times[j] >= o )
                        {
                            break;
                        }

                        // increase incomplete sampling psi
                        double dt = std::min(o, t_0) - std::max(std::max(min_age, d), times[j]);
                        psi_y_o += fossil[j]*dt;

                        size_t k = 0;
                        // increase running psi total for each observation
                        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                        {
                            if ( Fi->first.getMax() > times[j] )
                            {
                                double dt = 1.0;

                                // only compute dt if this is a non-singleton
                                if ( !( Fi->first.getMin() <= t_0 && (Fi->first.getMin() == Fi->first.getMax()) ) )
                                {
                                    dt = std::min(std::min(Fi->first.getMax(), o), t_0) - std::max(std::max(Fi->first.getMin(), d), times[j]);
                                }

                                psi[k] += fossil[j] * dt;
                            }
                        }
                    }

                    // include instantaneous sampling density
                    Psi[i] = log(fossil[oi]) + extended ? 0.0 : log(fossil[di]);

                    int count = 0;
                    double recip_o = 0.0;
                    double recip_y = 0.0;
                    double recip_oy = 0.0;

                    size_t k = 0;
                    // compute factors of the sum over each possible oldest/youngest observation
                    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                    {
                        count += Fi->second;

                        // compute sum of reciprocal oldest ranges
                        if ( Fi->first.getMax() >= o )
                        {
                            recip_o += Fi->second / psi[k];
                        }
                        // compute sum of reciprocal youngest ranges
                        if ( Fi->first.getMin() <= d )
                        {
                            recip_y += Fi->second / psi[k];

                            // compute sum of reciprocal oldest+youngest ranges
                            if ( Fi->first.getMax() >= o )
                            {
                                double f = Fi->second / psi[k];

                                recip_oy += f*(f-1.0);
                            }
                        }

                        // compute product of ranges
                        Psi[i] += log(psi[k]) * Fi->second;
                    }

                    // sum over each possible oldest/youngest observation
                    Psi[i] += extended ? log(recip_o) : log(recip_o * recip_y - recip_oy);

                    if ( complete == true )
                    {
                        // compute poisson density for count
                        Psi[i] -= RbMath::lnFactorial(count);
                    }
                    else
                    {
                        // compute poisson density for count + kappa, kappa >= 0
                        Psi[i] -= log(count);
                        Psi[i] += psi_y_o;

                        k = extended ? std::max(1, count-2) : count-1;

                        Psi[i] -= k*log(psi_y_o);
                        Psi[i] += k > 0 ? log(RbMath::incompleteGamma(psi_y_o, k, true, true)) : 0.0;
                    }
                }
                // only one fossil age
                else
                {
                    // include instantaneous sampling density
                    Psi[i] = ages.begin()->second * log(fossil[oi]);
                }
            }

            partial_likelihood[i] += Psi[i];
        }

        lnProbTimes += partial_likelihood[i];
    }

    // (the origin is not a speciation event)
    lnProbTimes -= log( birth[findIndex(origin)] );

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


/**
 * return the index i so that t_{i-1} > t >= t_i
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 is origin
 * t_l = 0.0
 */
size_t AbstractFossilizedBirthDeathProcess::findIndex(double t) const
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
    
    if ( t == 0.0 ) return 0.0;
    
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
void AbstractFossilizedBirthDeathProcess::redrawAge(size_t i, bool force)
{
    if ( force || GLOBAL_RNG->uniform01() < resampling )
    {
        dirty_taxa[i] = true;
        dirty_psi[i]  = true;

        age[i] = GLOBAL_RNG->uniform01()*(fbd_taxa[i].getMaxAge() - o_i[i]) + o_i[i];
    }
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
    age = stored_age;
    Psi = stored_Psi;

    dirty_psi  = std::vector<bool>(fbd_taxa.size(), false);
    dirty_taxa = std::vector<bool>(fbd_taxa.size(), false);

    touched = false;
}


void AbstractFossilizedBirthDeathProcess::touchSpecialization(DagNode *toucher, bool touchAll)
{
    if ( touched == false )
    {
        stored_likelihood = partial_likelihood;
        stored_age = age;
        stored_Psi = Psi;

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
void AbstractFossilizedBirthDeathProcess::prepareProbComputation()
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
        birth = heterogeneous_psi->getValue();
    }
    if ( homogeneous_psi != NULL )
    {
        fossil = std::vector<double>(num_intervals, homogeneous_psi->getValue() );
    }
    else
    {
        birth = heterogeneous_psi->getValue();
    }

    times = timeline->getValue();

    if ( times.size() < num_intervals )
    {
        times.insert(times.begin(), 0.0);
    }

    for (size_t interval = num_intervals; interval > 0; interval--)
    {
        // start with youngest interval
        size_t i = interval - 1;

        double ti = times[i];
        double b = birth[i];
        double d = death[i];
        double f = fossil[i];

        if (i > 0)
        {
            double r = (i == num_intervals - 1 ? homogeneous_rho->getValue() : 0.0);
            double t = times[i-1];

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
