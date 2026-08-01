#include "AbstractFossilizedBirthDeathRangeProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <sstream>
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
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 */
AbstractFossilizedBirthDeathRangeProcess::AbstractFossilizedBirthDeathRangeProcess(const DagNode *inspeciation,
                                                                         const DagNode *inextinction,
                                                                         const DagNode *inpsi,
                                                                         const TypedDagNode<double> *inrho,
                                                                         const TypedDagNode< RbVector<double> > *intimes,
                                                                         const std::string &incondition,
                                                                         const std::vector<Taxon> &intaxa,
                                                                         bool comp,
                                                                         const TypedDagNode<double> *inorigin,
                                                                         TypedDistribution<double> *inoriginprior) :
    taxa(intaxa),
    condition(incondition),
    homogeneous_rho(inrho),
    timeline( intimes ),
    origin_age( inorigin ),
    origin_prior( inoriginprior ),
    origin(0.0),
    max_birth(0),
    resampled(false),
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
    heterogeneous_lambda    = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda      = dynamic_cast<const TypedDagNode<double >*>(inspeciation);
    heterogeneous_mu        = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
    homogeneous_mu          = dynamic_cast<const TypedDagNode<double >*>(inextinction);
    heterogeneous_psi       = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inpsi);
    homogeneous_psi         = dynamic_cast<const TypedDagNode<double >*>(inpsi);

    // add the parameters to the model
    range_parameters.push_back( timeline );
    range_parameters.push_back( origin_age );

    // the prior's own parameters have to reach the DAG
    if ( origin_prior != NULL )
    {
        const std::vector<const DagNode*> &pars = origin_prior->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); it++)
        {
            range_parameters.push_back( *it );
        }
    }
    range_parameters.push_back( homogeneous_rho );
    range_parameters.push_back( homogeneous_lambda );
    range_parameters.push_back( heterogeneous_lambda );
    range_parameters.push_back( homogeneous_mu );
    range_parameters.push_back( heterogeneous_mu );
    range_parameters.push_back( homogeneous_psi );
    range_parameters.push_back( heterogeneous_psi );

    // setup the timeline
    if ( timeline == NULL )
    {
        num_intervals = 1;
    }
    else
    {
        num_intervals = timeline->getValue().size() + (timeline->getValue().front() != 0.0);
    }

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

    if( heterogeneous_lambda != NULL || heterogeneous_mu != NULL || heterogeneous_psi != NULL)
    {
        if ( timeline == NULL ) throw(no_timeline_err);

        if( heterogeneous_lambda != NULL )
        {
            num_rates = heterogeneous_lambda->getValue().size();
        }
        if( heterogeneous_mu != NULL )
        {
            if ( num_rates == 0 ) num_rates = heterogeneous_mu->getValue().size();

            if ( heterogeneous_mu->getValue().size() != num_rates ) throw(inconsistent_rates);
        }
        if( heterogeneous_psi != NULL )
        {
            if ( num_rates == 0 ) num_rates = heterogeneous_psi->getValue().size();

            if ( heterogeneous_psi->getValue().size() != num_rates ) throw(inconsistent_rates);
        }
    }
    else
    {
        num_rates = 1;
    }

    if ( num_rates != num_intervals )
    {
        std::stringstream ss;
        ss << "Number of rates does not match number of time intervals in fossilized birth death process.";
        throw(RbException(ss.str()));
    }

    ranges = std::vector<RangeEntry>(taxa.size());

    p_i         = std::vector<double>(num_intervals, 1.0);
    pS_i        = std::vector<double>(num_intervals, 1.0);
    q_i         = std::vector<double>(num_intervals, 0.0);
    q_tilde_i   = std::vector<double>(num_intervals, 0.0);

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    partial_likelihood = std::vector<double>(taxa.size(), 0.0);


    dirty_taxa  = std::vector<bool>(taxa.size(), true);

    record_complete = comp;

    updateRecord();

    prepareProbComputation();

    if ( times.front() > max_present_age )
    {
        throw(RbException("Timeline start time is older than youngest fossil occurrence."));
    }
}


/**
 * Derive everything the likelihood reads off the occurrence record: the per-taxon counts, which
 * reporting model each taxon falls under, the bounds first_min and last_max that constrain the augmented
 * extremes, and the default tau_K. Called from the constructor and again whenever a clamped
 * record replaces the occurrences.
 */
void AbstractFossilizedBirthDeathRangeProcess::updateRecord( void )
{
    // the record fixes the bounds and the reporting model; the sampled ages survive a re-derive
    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        ranges[i].first_min = 0.0;
        ranges[i].first_max = RbConstants::Double::inf;
        ranges[i].last_min  = 0.0;
        ranges[i].last_max  = RbConstants::Double::inf;
        ranges[i].singleton = true;
        ranges[i].record.clear();
    }

    max_present_age = RbConstants::Double::inf;

    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        ranges[i].record = taxa[i].getOccurrences();

        size_t count = 0;
        for ( size_t k = 0; k < ranges[i].record.size(); k++ )
        {
            const TimeInterval &bin = ranges[i].record[k].first;

            count += ranges[i].record[k].second;

            // find the oldest minimum age
            ranges[i].first_min = std::max(bin.getMin(), ranges[i].first_min);
            // find the youngest maximum age
            ranges[i].last_max = std::min(bin.getMax(), ranges[i].last_max);

            max_present_age = std::min(max_present_age, ranges[i].last_max);
        }

        ranges[i].singleton = ( count < 2 );

        ranges[i].first_max = taxa[i].getMaxAge();
        ranges[i].last_min  = taxa[i].getMinAge();

        // default the augmented youngest age to the youngest maximum (only resampled,
        // and only used, when the record has two extremes to order)
        ranges[i].last = ranges[i].last_max;
    }

}


/**
 * Adopt a clamped record's occurrences as the data. The augmented extremes and the value were drawn
 * against the old bins, so re-derive the record model and hand the value to the derived process to
 * be moved back into the new support.
 */
void AbstractFossilizedBirthDeathRangeProcess::setOccurrences( const std::vector<Taxon> &t )
{
    // clamping the record the process was built with is the ordinary case, and it must leave the
    // augmented ages the constructor drew against those occurrences exactly as they are
    bool changed = ( t.size() != taxa.size() );

    for (size_t i = 0; i < taxa.size() && changed == false; ++i)
    {
        if ( taxa[i].getOccurrences() != t[i].getOccurrences() ) changed = true;
        if ( taxa[i].isExtinct()     != t[i].isExtinct()     ) changed = true;
    }

    if ( changed == false ) return;

    taxa = t;

    updateRecord();

    prepareProbComputation();

    if ( times.front() > max_present_age )
    {
        throw(RbException("dnFossilRecord: the clamped record has an occurrence younger than the timeline start."));
    }

    dirty_taxa = std::vector<bool>(taxa.size(), true);

    // the extremes were drawn under the old reporting model, and which of them the new one even
    // instantiates may differ, so draw them again rather than patching the old values
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        drawAugmentedAges(i);
    }

    // the new record may still leave the value out of support, which no clipping repairs: the
    // range ends and births were drawn against the old bins. The MCMC redraws every free node
    // jointly when a start is -inf, so leave that to it rather than redrawing this one alone.
    clipAugmentedAges();
}


/**
 * Move any augmented extreme that the current ranges put outside its bin back inside it. A taxon
 * whose bin the ranges cannot accommodate at all is left alone, since only a different value fixes it.
 */
void AbstractFossilizedBirthDeathRangeProcess::clipAugmentedAges( void )
{
    prepareProbComputation();
    updateRanges();

    for (size_t i = 0; i < taxa.size(); i++)
    {
        ranges[i].clip();
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double AbstractFossilizedBirthDeathRangeProcess::computeLnProbabilityRanges( bool force )
{
    // prepare the probability computation
    prepareProbComputation();

    updateRanges();

    // the origin is the oldest birth; a supplied one pins it
    if ( origin_age != NULL && ranges[max_birth].birth != origin )
    {
        return RbConstants::Double::neginf;
    }

    // variable declarations and initialization
    double lnProb = 0.0;

    size_t num_rho_sampled = 0;
    size_t num_rho_unsampled = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        double b = ranges[i].birth;
        double d = ranges[i].death;

        double o = ranges[i].first;

        double max_age = taxa[i].getMaxAge();
        double min_age = taxa[i].getMinAge();

        double present = times.front();

        // check model constraints. tau_K sits between the range end and tau_1, which the
        // reporting term only checks for taxa whose two extremes are separate variables; a move
        // on the value can put it anywhere, so the support is enforced here for every taxon
        if ( ranges[i].isOrdered( present ) == false )
        {
            return RbConstants::Double::neginf;
        }
        // The status flag is the rho-sampling datum; the ranges are the psi-sampling data.
        // Seeing a taxon at the present pins its death there, but NOT seeing one leaves d
        // free: an unsampled lineage may still have survived, and pays 1-rho if it did.
        if ( taxa[i].isExtinct() == false && d != present )
        {
            return RbConstants::Double::neginf;
        }

        num_rho_sampled   += ( taxa[i].isExtinct() == false );          // l
        // a marginalized range closes with p(), which already carries the survived-unseen branch
        num_rho_unsampled += ( taxa[i].isExtinct() && d == present && marginalizesExtinction() == false );    // n - m - l

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
            for (size_t j = oi; j < bi; j++)
            {
                partial_likelihood[i] += q_i[j];
            }

            // skip the rest for extant taxa with no fossil samples
            if ( max_age == present )
            {
                lnProb += partial_likelihood[i];

                continue;
            }

            // replace q terms at oldest occurrence
            partial_likelihood[i] += q(oi, o, true) - q(oi, o);

            // include intermediate q_tilde terms
            for (size_t j = di; j < oi; j++)
            {
                partial_likelihood[i] += q_tilde_i[j];
            }

            // divide by q_tilde at the death time
            partial_likelihood[i] -= q( di, d, true);

            // close the range: an extinction density, or whatever the derived process puts there
            // when the extinction time is integrated out
            if ( d > present )
            {
                partial_likelihood[i] += rangeEndTerm( i, di, d );
            }

        }

        lnProb += partial_likelihood[i];
    }

    size_t ori = findIndex(origin);

    // the origin is not a speciation event
    lnProb -= log( birth[ori] );

    // a supplied prior applies to the oldest birth, which is the origin
    if ( origin_prior != NULL )
    {
        origin_prior->setValue( new double(origin) );
        lnProb += origin_prior->computeLnProbability();
    }

    // Extant tip age terms. Status is data
    double rho = homogeneous_rho->getValue();

    if ( num_rho_sampled > 0 )                          // seen at the present, so rho > 0
    {
        if ( rho == 0.0 ) return RbConstants::Double::neginf;
        lnProb += num_rho_sampled * log( rho );
    }
    if ( num_rho_unsampled > 0 )                        // survived unseen, so rho < 1
    {
        if ( rho == 1.0 ) return RbConstants::Double::neginf;
        lnProb += num_rho_unsampled * log( 1.0 - rho );
    }

    // condition on sampling
    if ( condition == "sampling" )
    {
        lnProb -= log( 1.0 - p(ori, origin, false) );
    }
    // condition on survival
    else if ( condition == "survival" )
    {
        lnProb -= log( 1.0 - p(ori, origin, true) );
    }

    if ( RbMath::isFinite(lnProb) == false )
    {
        return RbConstants::Double::neginf;
    }

    return lnProb;
}


// Total fossil-occurrence (reporting) log-density for a standalone reporting node
// (dnFossilRecord) conditioned on this range process. Self-contained: refreshes the piecewise-
// rate cache and per-taxon start/end times, then sums the per-taxon reporting terms --
double AbstractFossilizedBirthDeathRangeProcess::computeLnFossilTotal()
{
    prepareProbComputation();
    updateRanges();

    double lnProb = 0.0;
    for ( size_t i = 0; i < taxa.size(); ++i )
    {
        // an extant taxon with no fossil sample carries no reporting term; computeLnProbabilityRanges
        // skips it as well, so the fused and split forms omit the same taxa
        if ( taxa[i].getMaxAge() == times.front() ) continue;

        double r = computeLnFossilRecord(i);
        if ( r == RbConstants::Double::neginf )
        {
            return RbConstants::Double::neginf;
        }
        lnProb += r;
    }
    return lnProb;
}


// Fossil-occurrence (reporting) log-term for taxon i: the Psi[i] block factored out of
// computeLnProbabilityRanges (range/reporting split). Behavior-preserving -- returns
// exactly what was assigned to Psi[i] inline before the split.
double AbstractFossilizedBirthDeathRangeProcess::computeLnFossilRecord( size_t i ) const
{
    double d = ranges[i].death;
    double o = ranges[i].first;                             // oldest augmented age
    size_t oi = findIndex(o);
    double y  = ranges[i].last;                             // youngest augmented age
    size_t yi = findIndex(y);

    // with tau_K instantiated the interior spans (tau_K, tau_1); otherwise it runs down to the
    // range end, since an unreported occurrence may lie anywhere above it
    double min_age = ranges[i].singleton == false ? y : d;
    double max_age = o;

    double result = 0.0;

    const std::vector<std::pair<TimeInterval, size_t> > &ages = ranges[i].record;

    // if there is a range of fossil ages
    if ( min_age != max_age )
    {
        double Lambda = 0.0;                            // interior sampling rate over (last, first)

        std::vector<double> psi(ages.size(), 0.0);

        for (size_t j = 0; j < num_intervals; j++)
        {
            double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

            if ( t_0 <= min_age )
            {
                continue;
            }
            if ( times[j] >= max_age )
            {
                break;
            }

            // the kappa interior span
            double int_hi = std::min(max_age, t_0);
            double int_lo = std::max(min_age, times[j]);
            double dti = int_hi - int_lo;
            if ( dti > 0.0 ) Lambda += fossil[j]*dti;

            size_t k = 0;
            // increase running psi total for each observation
            for ( std::vector<std::pair<TimeInterval, size_t> >::const_iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
            {
                if ( Fi->first.getMin() < t_0 && Fi->first.getMax() > times[j] )
                {
                    double dt = 1.0;

                    // only compute dt if this is a non-singleton
                    if ( Fi->first.getMin() != Fi->first.getMax() )
                    {
                        // observed occurrence sampling rate over its interval
                        dt = std::min(Fi->first.getMax(), int_hi) - std::max(Fi->first.getMin(), int_lo);
                    }

                    psi[k] += fossil[j] * dt;
                }
            }
        }

         // instantaneous rate of the oldest specimen at tau1
        result = log(fossil[oi]);

        size_t count = 0;
        for ( size_t k = 0; k < ages.size(); k++ ) count += ages[k].second;
        double recip_old = 0.0;                         // sum_{i: tau1 in F_i} count_i / Psi(F_i)
        double recip_young = 0.0;                       // sum_{i: tauK in F_i} count_i / Psi(F_i)
        double diag = 0.0;

        size_t k = 0;
        // compute factors of the sum over each possible oldest/youngest observation
        for ( std::vector<std::pair<TimeInterval, size_t> >::const_iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
        {
            bool eligible_oldest = ( Fi->first.getMax() >= o );
            bool eligible_youngest = ( Fi->first.getMin() <= y );

            // compute sum of reciprocal oldest ranges
            if ( eligible_oldest )
            {
                recip_old += Fi->second / psi[k];
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
            result += log(psi[k]) * Fi->second;
        }

        // sum over which observation supplies the oldest specimen at tau1
        result += log(recip_old) - RbMath::lnFactorial(int(count));

        if ( ranges[i].singleton == false )
        {
            if ( !( o >= y && y >= d && y <= ranges[i].last_max ) )
            {
                return RbConstants::Double::neginf;
            }
            // youngest instantaneous density + sum over which observation is the youngest,
            // excluding the diagonal where a single occurrence supplies both extremes
            result += log(fossil[yi]) + log(recip_young - diag/recip_old);

            // the first/last rule leaves the interior occurrences unreported, so their count
            // is marginalized; a complete record has none to marginalize
            if ( record_complete == false )
            {
                double S1 = 0.0, f = 1.0;
                for ( size_t kap = 0; kap < 200; kap++ )
                {
                    S1 += f;
                    f *= Lambda / double(count - 1 + kap);
                    if ( f < 1e-16 * S1 ) break;
                }
                result += log(S1);
            }
        }
    }
    // only one fossil age
    else
    {
        // include instantaneous sampling density
        result = ages.begin()->second * log(fossil[oi]);
    }

    return result;
}


/**
 * return the index i so that t_{i-1} > t >= t_i
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 is origin
 * t_l = 0.0
 */
size_t AbstractFossilizedBirthDeathRangeProcess::findIndex(double t) const
{
    return std::prev(std::upper_bound( times.begin(), times.end(), t)) - times.begin();
}


/**
 * p_i(t)
 */
double AbstractFossilizedBirthDeathRangeProcess::p( size_t i, double t, bool survival ) const
{
    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = survival ? 0.0 : fossil[i];
    double pi = survival ? pS_i[i] : p_i[i];
    double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];
    
    double diff = b - d - f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*b*f);

    // survival takes f = 0, so b == d leaves A = 0, where the closed form is singular
    if ( A < 1E-10 )
    {
        double s = 1.0 - (1.0 - r)*pi;

        return 1.0 - s/(1.0 + b*s*dt);
    }

    double B = ( (1.0 - 2.0*(1.0-r)*pi )*b + d + f ) / A;

    double ln_e = -A*dt;

    double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

    return (b + d + f - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);
}


/**
 * q_i(t)
 */
double AbstractFossilizedBirthDeathRangeProcess::q( size_t i, double t, bool tilde ) const
{
    if ( t == times[i] ) return 0.0;
    
    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
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
/**
 * The augmented first (tau_1) and last (tau_K) ages, which no monitor can otherwise reach.
 * A record with one occurrence or an unreported youngest has last == first.
 */
void AbstractFossilizedBirthDeathRangeProcess::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    if ( n == "getAugmentedFirstAges" || n == "getAugmentedLastAges" || n == "getBirthAges" || n == "getDeathAges" )
    {
        // range_start/range_end are a byproduct of the density, so refresh them from the value
        AbstractFossilizedBirthDeathRangeProcess *self = const_cast<AbstractFossilizedBirthDeathRangeProcess *>( this );
        self->updateRanges();

        // with the extinction times marginalized out there is none to report: range_end holds the range
        // end instead, and a sampled ancestor's is jointly distributed with the unobserved
        // speciation separating it from its descendant, so it is not recoverable here
        if ( n == "getDeathAges" && marginalizesExtinction() == true )
        {
            throw RbException("getDeathAges is unavailable when extended=false: the extinction times are marginalized out rather than sampled. Use getAugmentedLastAges for the youngest occurrence ages, or extended=true to sample extinction times.");
        }

        rv.clear();
        for (size_t i = 0; i < ranges.size(); i++)
        {
            rv.push_back( n == "getAugmentedFirstAges" ? ranges[i].first :
                        ( n == "getAugmentedLastAges"  ? ranges[i].last  :
                        ( n == "getBirthAges"          ? ranges[i].birth : ranges[i].death ) ) );
        }
    }
    else
    {
        throw RbException() << "The fossilized birth death range process does not have a member method called '" << n << "'.";
    }
}


/**
 * The origin of the process, which is the oldest birth.
 */
void AbstractFossilizedBirthDeathRangeProcess::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, double &rv) const
{
    if ( n == "getOrigin" )
    {
        // origin is a byproduct of the density and is not restored, so refresh it from the value
        AbstractFossilizedBirthDeathRangeProcess *self = const_cast<AbstractFossilizedBirthDeathRangeProcess *>( this );
        self->updateRanges();

        rv = origin;
    }
    else
    {
        throw RbException() << "The fossilized birth death range process does not have a member method called '" << n << "'.";
    }
}


/**
 * Draw the augmented extremes for taxon i, nested (range_end <= last <= first), under the current
 * reporting model.
 */
void AbstractFossilizedBirthDeathRangeProcess::drawAugmentedAges(size_t i)
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    bool augment_youngest = ranges[i].singleton == false;

    double lo, hi;

    // youngest augmented age, at or above the death and within its reported bin. A non-extended
    // range ends at this age rather than at an extinction, so the death does not bound it.
    if ( augment_youngest )
    {
        lo = marginalizesExtinction() ? taxa[i].getMinAge() : std::max( ranges[i].death, taxa[i].getMinAge() );
        hi = ranges[i].last_max;
        ranges[i].last = ( hi > lo ) ? rng->uniform01()*(hi - lo) + lo : lo;
    }
    else
    {
        ranges[i].last = ranges[i].death;
    }

    // oldest augmented age, at or above the youngest and the oldest reported minimum
    lo = std::max( ranges[i].last, ranges[i].first_min );
    hi = ranges[i].first_max;
    ranges[i].first = ( hi > lo ) ? rng->uniform01()*(hi - lo) + lo : lo;

    // a single occurrence is its own youngest
    if ( augment_youngest == false ) ranges[i].last = ranges[i].first;

    // a non-extended range ends at its youngest augmented age, so the initial tree is built there.
    // An extant range still ends at the present, and its tau_K is a fossil age, not its tip.
    if ( marginalizesExtinction() == true && taxa[i].isExtinct() == true ) ranges[i].death = ranges[i].last;

}

// The augmented ages move only through mvResampleAugmentedAges. Without it they stay at their
// initial draw and the chain silently samples the wrong space, so say so once at startup.
void AbstractFossilizedBirthDeathRangeProcess::warnIfNoResampleMove( void ) const
{
    if ( has_resample_move == false && warned_no_resample == false )
    {
        warned_no_resample = true;
        RBOUT("Warning: no mvResampleAugmentedAges move; augmented ages will not be sampled.");
    }
}


/**
 * This process omits every fossil-occurrence density, so on its own it
 * is not a density over the record at all and psi is left with no data. A dnFossilRecord supplies
 * that term; warn once if none does.
 */
void AbstractFossilizedBirthDeathRangeProcess::warnIfNoReportingNode( void ) const
{
    if ( has_reporting_node == false && warned_no_reporting == false )
    {
        warned_no_reporting = true;
        RBOUT("Warning: no dnFossilRecord node; fossil sampling is not scored and psi has no data.");
    }
}


void AbstractFossilizedBirthDeathRangeProcess::resampleFirstLast(size_t i)
{
    stored_first.resize( taxa.size() );
    stored_last.resize( taxa.size() );
    for (size_t k = 0; k < taxa.size(); ++k)
    {
        stored_first[k] = ranges[k].first;
        stored_last[k]  = ranges[k].last;
    }
    resampled = true;

    // a non-extended extinct range ends at its tip, which a move samples and updateRanges
    // reads back. With one occurrence that tip is tau_1 too, so there is nothing to draw at all
    if ( marginalizesExtinction() == true && taxa[i].isExtinct() == true && ranges[i].singleton )
    {
        return;
    }

    // the bin that reported it, and nothing state-dependent: this is an independence proposal
    // whose Hastings ratio is 1 only while its support is fixed by the data
    double hi = ranges[i].first_max;
    ranges[i].first = GLOBAL_RNG->uniform01()*(ranges[i].first_min - hi) + hi;

    // an extinct non-extended tip is the augmented youngest age itself, so tau_K is not drawn
    // here. An extant tip sits at the present instead, so its tau_K still is.
    if ( marginalizesExtinction() == true && taxa[i].isExtinct() == true )
    {
        return;
    }

    // a single occurrence is its own youngest, as is an exchangeable record whose true youngest
    // may be unreported
    if ( ranges[i].singleton == false )
    {
        ranges[i].last = GLOBAL_RNG->uniform01()*(ranges[i].last_max - taxa[i].getMinAge()) + taxa[i].getMinAge();
    }
    else
    {
        ranges[i].last = ranges[i].first;
    }

}


void AbstractFossilizedBirthDeathRangeProcess::drawRanges()
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // bracket birth times: an unbounded max_age would push every birth to infinity, so fall
    // back to the oldest lower bound the taxon's occurrences provide
    double max = 0;
    for (size_t i = 0; i < taxa.size(); i++)
    {
        double o = taxa[i].getMaxAge();
        if ( RbMath::isFinite(o) == false )
        {
            o = 0.0;
            for ( size_t k = 0; k < ranges[i].record.size(); k++ )
            {
                o = std::max( o, ranges[i].record[k].first.getMin() );
            }
        }
        if ( o > max ) max = o;
    }
    max *= 1.1;
    if ( max == 0.0 ) max = 1.0;

    double present = times.front();

    for (size_t i = 0; i < taxa.size(); i++)
    {
        // Draw d over its range, then the augmented ages NESTED (d <= last <= first)
        ranges[i].death = taxa[i].isExtinct() ? rng->uniform01()*(ranges[i].last_max - present) + present : present;
        ranges[i].birth = max;

        drawAugmentedAges(i);
    }

    // place births oldest-first, each inside a lineage already placed and still alive at it
    std::vector<size_t> order( taxa.size() );
    for (size_t i = 0; i < taxa.size(); i++) order[i] = i;

    std::sort( order.begin(), order.end(), [&](size_t a, size_t b) { return ranges[a].first > ranges[b].first; } );

    for (size_t k = 0; k < order.size(); k++)
    {
        size_t i = order[k];

        // the oldest birth is the origin: a supplied one pins it, a prior supplies its support
        if ( k == 0 )
        {
            if ( origin_age != NULL )
            {
                ranges[i].birth = origin_age->getValue();
            }
            else if ( origin_prior != NULL )
            {
                // the origin has to clear every oldest age, which a blind draw rarely does
                double oldest = ranges[i].first;

                origin_prior->redrawValue();
                for (size_t t = 0; t < 1000 && origin_prior->getValue() <= oldest; t++)
                {
                    origin_prior->redrawValue();
                }

                ranges[i].birth = origin_prior->getValue();
            }
            else
            {
                ranges[i].birth = ranges[i].first + rng->uniform01()*(max - ranges[i].first);
            }

            continue;
        }

        size_t pick = size_t( rng->uniform01()*k );
        if ( pick >= k ) pick = k - 1;

        size_t a = order[pick];

        // sorted oldest-first, so the window is never empty
        double lo = std::max( ranges[i].first, ranges[a].death );

        ranges[i].birth = lo + rng->uniform01()*(ranges[a].birth - lo);
    }
}


void AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(const DagNode *toucher)
{
    dirty_taxa = std::vector<bool>(taxa.size(), false);

    resampled = false;
    touched = false;
}


void AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(const DagNode *toucher)
{
    partial_likelihood = stored_likelihood;

    if ( resampled )
    {
        for (size_t i = 0; i < taxa.size(); ++i)
        {
            ranges[i].first = stored_first[i];
            ranges[i].last  = stored_last[i];
        }
    }

    dirty_taxa = std::vector<bool>(taxa.size(), false);

    resampled = false;
    touched = false;
}


void AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    if ( touched == false )
    {
        stored_likelihood = partial_likelihood;

        dirty_taxa = std::vector<bool>(taxa.size(), true);

        if ( toucher == timeline || toucher == homogeneous_psi || toucher == heterogeneous_psi || touchAll )
        {
        }
    }

    touched = true;
}


/**
 *
 *
 */
void AbstractFossilizedBirthDeathRangeProcess::prepareProbComputation( void ) const
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

    for (size_t i = 0; i < num_intervals; i++)
    {
        double ti = times[i];
        double b = birth[i];
        double d = death[i];
        double f = fossil[i];

        if (i < num_intervals-1)
        {
            double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
            double t = times[i+1];

            double diff = b - d - f;
            double dt   = t - ti;

            double A = sqrt( diff*diff + 4.0*b*f);
            double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

            double ln_e = -A*dt;

            double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

            q_i[i]       = log(4.0) + ln_e - 2.0*log(tmp);
            q_tilde_i[i] = 0.5 * ( q_i[i] - (b+d+f)*dt );
            p_i[i+1]       = (b + d + f - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);

            if ( condition == "survival" )
            {
                // survival ignores fossil sampling, so this is the f = 0 case of the recursion
                // above, carried on its own probability rather than the sampled one
                diff = b - d;

                A = fabs( diff );

                double s = 1.0 - (1.0 - r)*pS_i[i];

                if ( A < 1E-10 )
                {
                    // b == d leaves A = 0, where the closed form is singular but its limit is not
                    pS_i[i+1] = 1.0 - s/(1.0 + b*s*dt);
                }
                else
                {
                    B = ( (1.0 - 2.0*(1.0-r)*pS_i[i] )*b + d ) / A;

                    ln_e = -A*dt;

                    tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

                    pS_i[i+1] = (b + d - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);
                }
            }
        }
    }
}


/**
 * Swap the parameters held by this distribution.
 * 
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
    else if (oldP == origin_age)
    {
        origin_age = static_cast<const TypedDagNode<double>* >( newP );
    }
}
