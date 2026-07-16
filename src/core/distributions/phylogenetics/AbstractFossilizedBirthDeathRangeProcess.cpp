#include "AbstractFossilizedBirthDeathRangeProcess.h"

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
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 * \param[in]    re             Augmented age resampling weight.
 */
AbstractFossilizedBirthDeathRangeProcess::AbstractFossilizedBirthDeathRangeProcess(const DagNode *inspeciation,
                                                                         const DagNode *inextinction,
                                                                         const DagNode *inpsi,
                                                                         const TypedDagNode<double> *inrho,
                                                                         const TypedDagNode< RbVector<double> > *intimes,
                                                                         const std::string &incondition,
                                                                         const std::vector<Taxon> &intaxa,
                                                                         const std::string &s,
                                                                         bool re,
                                                                         const TypedDagNode<double> *inorigin) :
    taxa(intaxa),
    condition(incondition),
    homogeneous_rho(inrho),
    timeline( intimes ),
    origin_age( inorigin ),
    origin(0.0),
    reporting(s),
    resampled(false),
    resampling(re),
    touched(false),
    report_internally(true)
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

    b_i = std::vector<double>(taxa.size(), 0.0);
    d_i = std::vector<double>(taxa.size(), 0.0);
    o_i = std::vector<double>(taxa.size(), 0.0);
    y_i = std::vector<double>(taxa.size(), RbConstants::Double::inf);

    p_i         = std::vector<double>(num_intervals, 1.0);
    pS_i        = std::vector<double>(num_intervals, 1.0);
    q_i         = std::vector<double>(num_intervals, 0.0);
    q_tilde_i   = std::vector<double>(num_intervals, 0.0);

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    partial_likelihood = std::vector<double>(taxa.size(), 0.0);

    first         = std::vector<double>(taxa.size(), 0.0);
    last   = std::vector<double>(taxa.size(), 0.0);
    Psi         = std::vector<double>(taxa.size(), 0.0 );

    dirty_taxa  = std::vector<bool>(taxa.size(), true);
    dirty_psi   = std::vector<bool>(taxa.size(), true);

    double max_present = RbConstants::Double::inf;

    max_count = 0;
    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        std::map<TimeInterval, size_t> ages = taxa[i].getOccurrences();
        size_t count = 0;
        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++ )
        {
            // find the oldest minimum age
            o_i[i] = std::max(Fi->first.getMin(), o_i[i]);
            // find the youngest maximum age
            y_i[i] = std::min(Fi->first.getMax(), y_i[i]);

            max_present = std::min(max_present, y_i[i]);
            count += Fi->second;
        }
        // the largest per-taxon occurrence count is the implicit reporting cap K
        // for the uniform model (see effectiveReporting)
        max_count = std::max(max_count, count);
        // default the augmented youngest age to the youngest maximum (only resampled,
        // and only used, under first/last conditioning)
        last[i] = y_i[i];
    }

    prepareProbComputation();

    if ( times.front() > max_present )
    {
        throw(RbException("Timeline start time is older than youngest fossil first."));
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

    updateStartEndTimes();

    // a supplied origin must be at least as old as every sampled birth
    if ( origin_age != NULL )
    {
        for (size_t i = 0; i < taxa.size(); ++i)
        {
            if ( origin < b_i[i] )
            {
                return RbConstants::Double::neginf;
            }
        }
    }

    // variable declarations and initialization
    double lnProb = 0.0;

    size_t num_rho_sampled = 0;
    size_t num_rho_unsampled = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        double b = b_i[i];
        double d = d_i[i];

        double o = first[i];

        double max_age = taxa[i].getMaxAge();
        double min_age = taxa[i].getMinAge();

        double present = times.front();

        // check model constraints
        //if ( !( b > max_age && min_age >= d && d >= 0.0 ) )
        if ( !( b > o && o >= d && o >= o_i[i] && y_i[i] >= d && d >= present ) )
        {
            return RbConstants::Double::neginf;
        }
        if ( (d > present) != taxa[i].isExtinct() )
        {
            return RbConstants::Double::neginf;
        }

        // count the number of rho-sampled tips
        num_rho_sampled   += (d == present && min_age == present);  // l
        num_rho_unsampled += (d == present && min_age > present); // n - m - l

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

            // include extinction density
            if ( d > present ) partial_likelihood[i] += log( death[di] );

            if ( report_internally )
            {
                if ( dirty_psi[i] || force )
                {
                    Psi[i] = computeLnFossilRecord(i);
                    if ( Psi[i] == RbConstants::Double::neginf )
                    {
                        return RbConstants::Double::neginf;
                    }
                }

                partial_likelihood[i] += Psi[i];
            }

            // Jacobian for the auto-resampled tau_1 ~ Uniform(lo, hi): under
            // u = (tau_1 - lo)/(hi - lo) the resample is symmetric on [0,1], and
            // log(hi - lo) is the change of variables. It is the only channel to the
            // acceptance ratio, the auto-resample being a touch side-effect with no
            // Hastings. Omitted when the augmentation is frozen (resample=false).
            if ( resampling == true )
            {
                double lo, hi;
                firstLastSupport(i, lo, hi);
                if ( hi > lo ) partial_likelihood[i] += log( hi - lo );
            }
        }

        lnProb += partial_likelihood[i];
    }

    size_t ori = findIndex(origin);

    // when the origin is not supplied, the oldest sampled birth is the process
    // origin and is not a speciation event
    if ( origin_age == NULL )
    {
        lnProb -= log( birth[ori] );
    }
    else
    {
        double max_birth = 0.0;
        for (size_t i = 0; i < taxa.size(); ++i)
        {
            max_birth = std::max(max_birth, b_i[i]);
        }

        size_t mbi = findIndex(max_birth);

        lnProb += q(ori, origin) - q(mbi, max_birth);

        for (size_t j = mbi; j < ori; ++j)
        {
            lnProb += q_i[j];
        }
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
// (dnFossilRecord) conditioned on this skeleton. Self-contained: refreshes the piecewise-
// rate cache and per-taxon start/end times, then sums the per-taxon reporting terms --
// exactly what computeLnProbabilityRanges adds inline when report_internally is true.
double AbstractFossilizedBirthDeathRangeProcess::computeLnFossilTotal()
{
    prepareProbComputation();
    updateStartEndTimes();

    double lnProb = 0.0;
    for ( size_t i = 0; i < taxa.size(); ++i )
    {
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
// computeLnProbabilityRanges (skeleton/reporting split). Behavior-preserving -- returns
// exactly what was assigned to Psi[i] inline before the split.
double AbstractFossilizedBirthDeathRangeProcess::computeLnFossilRecord( size_t i ) const
{
    double d = d_i[i];
    double o = first[i];
    double min_age = taxa[i].getMinAge();
    double max_age = taxa[i].getMaxAge();
    size_t oi = findIndex(o);

    double result = 0.0;

                std::map<TimeInterval, size_t> ages = taxa[i].getOccurrences();

                // if there is a range of fossil ages
                if ( min_age != max_age )
                {
                    if ( effectiveReporting(i) == "firstlast" )
                    {

                    double psi_int = 0.0;                            // interior sampling rate over (last, first)

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
                                    // observed occurrence sampling rate over its interval (below the oldest);
                                    // the youngest only enters via the interior rate psi_int, not here
                                    dt = std::min(std::min(Fi->first.getMax(), o), t_0) - std::max(std::max(Fi->first.getMin(), d), times[j]);
                                }

                                psi[k] += fossil[j] * dt;
                            }
                        }
                    }

                    // include instantaneous sampling density (oldest; the youngest is added
                    // inside the count >= 2 branch below, so a single-occurrence taxon - which
                    // skips that branch - never double-counts an instantaneous density)
                    result = log(fossil[oi]);

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
                        result += log(psi[k]) * Fi->second;
                    }

                    // sum over each possible oldest observation
                    result += log(recip);

                    if ( count >= 2 )
                    {
                        // Condition on the oldest and youngest occurrences: sum over which
                        // observation is the youngest (recip_young) and marginalize the
                        // kappa >= 0 unobserved interior specimens in (last, first). Inside
                        // the count >= 2 branch, so single-occurrence taxa ignore 'last'.
                        if ( !( o >= last[i] && last[i] >= d && last[i] <= y_i[i] && last[i] >= min_age ) )
                        {
                            return RbConstants::Double::neginf;
                        }
                        // youngest instantaneous density + sum over which observation is the youngest,
                        // excluding the diagonal where a single occurrence supplies both extremes
                        result += log(fossil[yi]) + log(recip_young - diag/recip);

                        double S1 = 0.0, f = 1.0;
                        for ( size_t kap = 0; kap < 200; kap++ )
                        {
                            S1 += f;
                            f *= psi_int / double(count - 1 + kap);
                            if ( f < 1e-16 * S1 ) break;
                        }
                        result += log(S1) - RbMath::lnFactorial(count);
                    }
                    // count == 1 (single occurrence, first == last): no interior term
                                    }
                    else  // "uniform" (exchangeable) or "complete"
                    {
                    // Truly-exchangeable incomplete sampling: the reported occurrences are a
                    // uniform subset, so the true oldest (tau1 = first[i]) may be unobserved and
                    // is augmented up to the birth. Unobserved specimens fall anywhere in
                    // (d, tau1), so Lambda integrates psi over that range, not the observed span.
                    double Lambda = 0.0;                            // Psi(d, tau1): sampling over the whole sampled interval
                    std::vector<double> psi(ages.size(), 0.0);      // Psi(F_i): per-occurrence integral over its interval, capped at [d,tau1]

                    for (size_t j = 0; j < num_intervals; j++)
                    {
                        double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

                        if ( t_0 <= d ) continue;
                        if ( times[j] >= o ) break;

                        double dL = std::min(o, t_0) - std::max(d, times[j]);
                        if ( dL > 0.0 ) Lambda += fossil[j]*dL;

                        size_t k = 0;
                        for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                        {
                            if ( Fi->first.getMin() < t_0 && Fi->first.getMax() > times[j] )
                            {
                                double dt = std::min(std::min(Fi->first.getMax(), o), t_0) - std::max(std::max(Fi->first.getMin(), d), times[j]);
                                if ( dt > 0.0 ) psi[k] += fossil[j] * dt;
                            }
                        }
                    }

                    // instantaneous rate of the oldest specimen at tau1
                    result = log(fossil[oi]);

                    int count = 0;
                    double recip_old = 0.0;                         // sum_{i: tau1 in F_i} count_i / Psi(F_i)

                    size_t k = 0;
                    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
                    {
                        count += Fi->second;
                        // occurrence i can be the labeled oldest specimen iff its interval contains tau1
                        if ( Fi->first.getMin() <= o && Fi->first.getMax() >= o )
                        {
                            recip_old += Fi->second / psi[k];
                        }
                        result += log(psi[k]) * Fi->second;         // log prod_i Psi(F_i)^{count_i}
                    }

                    if ( effectiveReporting(i) == "complete" )
                    {
                        result -= RbMath::lnFactorial(count);
                    }
                    else
                    {
                        // P(N >= count) and P(N >= count+1) for N ~ Poisson(Lambda), by upper-tail
                        // summation (no catastrophic cancellation).
                        double Pk = 0.0;
                        double pmf_count = exp( -Lambda + count*log(Lambda) - RbMath::lnFactorial(count) );
                        double t = pmf_count;
                        for ( int n = count; n < count + 100000; n++ )
                        {
                            Pk += t;
                            t *= Lambda / double(n+1);
                            if ( t < 1e-17 * Pk && n > (int)Lambda ) break;
                        }
                        double Pk1 = Pk - pmf_count;                // P(N >= count+1)

                        // bracket = P>=k * recip_old       (oldest specimen is one of the reported)
                        //         + ( P>=k - (k/Lambda) P>=k+1 )   (oldest unobserved; all reported interior)
                        double bracket = Pk*recip_old + Pk - (double(count)/Lambda)*Pk1;
                        if ( bracket <= 0.0 ) return RbConstants::Double::neginf;

                        // +Lambda: the q_tilde terms already carry e^{-Lambda}, which the
                        // tails re-introduce.
                        result += RbMath::lnFactorial(count) - count*log(Lambda) + log(bracket) + Lambda;
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
std::vector<double>& AbstractFossilizedBirthDeathRangeProcess::getAges(void)
{
    return first;
}


// Per-taxon reporting model. "uniform" reports a random subset capped at K = max_count:
// a taxon at the cap may have unreported fossils and gets the exchangeable
// marginalization, while one below it kept its whole record, which is the complete case.
std::string AbstractFossilizedBirthDeathRangeProcess::effectiveReporting(size_t i) const
{
    if ( reporting != "uniform" )
    {
        return reporting;
    }

    size_t count = 0;
    std::map<TimeInterval, size_t> ages = taxa[i].getOccurrences();
    for ( std::map<TimeInterval, size_t>::iterator Fi = ages.begin(); Fi != ages.end(); Fi++ )
    {
        count += Fi->second;
    }

    return ( count == max_count ) ? "uniform" : "complete";
}


// Set the reporting model. dnFossilRecord pushes its `reporting=` arg onto the skeleton via
// this setter, so the single `reporting` member is the one source of truth -- driving both the
// tau1 support (firstLastSupport) and the per-taxon reporting term (effectiveReporting).
void AbstractFossilizedBirthDeathRangeProcess::setReportingModel( const std::string &s )
{
    reporting = s;
}


// Support of the augmented oldest age tau_1, in one place so the resampling proposal
// and its Jacobian (computeLnProbabilityRanges) agree.
//  - "uniform": the true oldest may be unobserved and older than every reported
//    occurrence, so tau_1 ranges up to the birth b -- an hi that depends on b, which
//    is what the Jacobian corrects for.
//  - "firstlast"/"complete": the oldest observed occurrence IS the oldest fossil, so
//    tau_1 is bounded by its bin [o_i, max_age] and the Jacobian is a constant.
//  - an unbounded oldest occurrence (max_age = Inf) leaves tau_1 with no upper bound from the
//    data, so the process supplies the one it implies: no fossil can predate the birth.
void AbstractFossilizedBirthDeathRangeProcess::firstLastSupport(size_t i, double &lo, double &hi) const
{
    lo = std::max(o_i[i], d_i[i]);
    if ( effectiveReporting(i) == "uniform" )
    {
        hi = ( b_i[i] > lo ) ? b_i[i] : std::max(taxa[i].getMaxAge(), b_i[i]);
    }
    else
    {
        hi = std::max(taxa[i].getMaxAge(), lo);
    }

    if ( RbMath::isFinite(hi) == false )
    {
        hi = std::max(b_i[i], lo);
    }
}


void AbstractFossilizedBirthDeathRangeProcess::resampleFirstLast(size_t i)
{
    stored_first = first;
    stored_last = last;
    resampled = true;

    double _lo, _hi;
    firstLastSupport(i, _lo, _hi);
    first[i] = GLOBAL_RNG->uniform01()*(_hi - _lo) + _lo;

    // also augment the youngest occurrence age (single-occurrence taxa skip the
    // count >= 2 likelihood branch, so 'last' is inert there and the value drawn
    // here is simply unused)
    last[i] = GLOBAL_RNG->uniform01()*(y_i[i] - taxa[i].getMinAge()) + taxa[i].getMinAge();
}


void AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(const DagNode *toucher)
{
    dirty_psi  = std::vector<bool>(taxa.size(), false);
    dirty_taxa = std::vector<bool>(taxa.size(), false);

    resampled = false;
    touched = false;
}


void AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(const DagNode *toucher)
{
    partial_likelihood = stored_likelihood;
    Psi = stored_Psi;

    if ( resampled )
    {
        first = stored_first;
        last = stored_last;
    }

    dirty_psi  = std::vector<bool>(taxa.size(), false);
    dirty_taxa = std::vector<bool>(taxa.size(), false);

    resampled = false;
    touched = false;
}


void AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    if ( touched == false )
    {
        stored_likelihood = partial_likelihood;
        stored_Psi = Psi;

        dirty_taxa = std::vector<bool>(taxa.size(), true);

        if ( toucher == timeline || toucher == homogeneous_psi || toucher == heterogeneous_psi || touchAll )
        {
            dirty_psi  = std::vector<bool>(taxa.size(), true);
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
                diff = b - d;

                A = sqrt( diff*diff);
                B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d ) / A;

                ln_e = -A*dt;

                tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

                pS_i[i]  = (b + d - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);
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
