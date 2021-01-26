#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

#include "DistributionExponential.h"
#include "AbstractPiecewiseConstantFossilizedRangeProcess.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "DagNode.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
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
AbstractPiecewiseConstantFossilizedRangeProcess::AbstractPiecewiseConstantFossilizedRangeProcess(const DagNode *inspeciation,
                                                                                                 const DagNode *inextinction,
                                                                                                 const DagNode *inpsi,
                                                                                                 const DagNode *incounts,
                                                                                                 const TypedDagNode<double> *inrho,
                                                                                                 const TypedDagNode< RbVector<double> > *intimes,
                                                                                                 const std::vector<Taxon> &intaxa,
                                                                                                 bool ua) :
    ascending(false), homogeneous_rho(inrho), timeline( intimes ), fbd_taxa(intaxa), ages_from_counts(ua)
{
    // initialize all the pointers to NULL
    homogeneous_lambda             = NULL;
    homogeneous_mu                 = NULL;
    homogeneous_psi                = NULL;
    heterogeneous_lambda           = NULL;
    heterogeneous_mu               = NULL;
    heterogeneous_psi              = NULL;
    fossil_count                   = NULL;
    fossil_count_data              = NULL;

    RbException no_timeline_err = RbException("No time intervals provided for piecewise constant fossilized birth death process");

    heterogeneous_lambda = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda = dynamic_cast<const TypedDagNode<double >*>(inspeciation);

    range_parameters.push_back( homogeneous_lambda );
    range_parameters.push_back( heterogeneous_lambda );

    if ( heterogeneous_lambda == NULL && homogeneous_lambda == NULL)
    {
        throw(RbException("Speciation rate must be of type RealPos or RealPos[]"));
    }
    else if( heterogeneous_lambda != NULL )
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

    if ( heterogeneous_mu == NULL && homogeneous_mu == NULL)
    {
        throw(RbException("Extinction rate must be of type RealPos or RealPos[]"));
    }
    else if( heterogeneous_mu != NULL )
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

    if ( heterogeneous_psi == NULL && homogeneous_psi == NULL)
    {
        throw(RbException("Fossilization rate must be of type RealPos or RealPos[]"));
    }
    else if( heterogeneous_psi != NULL )
    {
        if( timeline == NULL ) throw(no_timeline_err);

        if (heterogeneous_psi->getValue().size() != timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of fossilization rates (" << heterogeneous_psi->getValue().size() << ") does not match number of time intervals (" << timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
    }

    fossil_count      = dynamic_cast<const TypedDagNode<long> *>(incounts);
    fossil_count_data = dynamic_cast<const TypedDagNode<HomologousDiscreteCharacterData<NaturalNumbersState> >*>(incounts);

    range_parameters.push_back( fossil_count );
    range_parameters.push_back( fossil_count_data );

    if ( fossil_count != NULL && homogeneous_psi == NULL)
    {
        throw(RbException("Heterogeneous fossil sampling rates provided, but homogeneous fossil counts"));
    }
    else if ( fossil_count_data != NULL )
    {
        if( timeline == NULL ) throw(no_timeline_err);

        if ( fossil_count_data->getValue().getNumberOfTaxa() != fbd_taxa.size())
        {
            std::stringstream ss;
            ss << "Number of species fossil counts (" << fossil_count_data->getValue().getNumberOfTaxa() << ") does not match number of taxa (" << fbd_taxa.size() << ")";
            throw(RbException(ss.str()));
        }
        else if ( fossil_count_data->getValue().getNumberOfCharacters() != timeline->getValue().size() + 1)
        {
            std::stringstream ss;
            ss << "Number of fossil counts per species (" << fossil_count_data->getValue().getNumberOfCharacters() << ") does not match number of time intervals (" << timeline->getValue().size() + 1 << ")";
            throw(RbException(ss.str()));
        }
    }

    range_parameters.push_back( homogeneous_rho );
    range_parameters.push_back( timeline );
    
    num_intervals = timeline == NULL ? 1 : timeline->getValue().size()+1;

    if ( num_intervals > 1 )
    {
        std::vector<double> times = timeline->getValue();
        std::vector<double> times_sorted_ascending = times;
        std::vector<double> times_sorted_descending = times;

        sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );
        sort(times_sorted_descending.rbegin(), times_sorted_descending.rend() );

        if( times == times_sorted_ascending )
        {
            ascending = true;
        }
        else if ( times != times_sorted_descending )
        {
            throw(RbException("Interval times must be provided in order"));
        }
    }

    b_i = std::vector<double>(fbd_taxa.size(), 0.0);
    d_i = std::vector<double>(fbd_taxa.size(), 0.0);

    H = std::vector<double>(fbd_taxa.size(), 0.0);

    p_i         = std::vector<double>(num_intervals, 1.0);
    q_i         = std::vector<double>(num_intervals, 0.0);
    q_tilde_i   = std::vector<double>(num_intervals, 0.0);

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    oldest_intervals = std::vector<size_t>( fbd_taxa.size(), 0 );
    youngest_intervals = std::vector<size_t>( fbd_taxa.size(), num_intervals - 1 );

    updateIntervals();
}

/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double AbstractPiecewiseConstantFossilizedRangeProcess::computeLnProbabilityRanges( void ) const
{
    // prepare the probability computation
    updateIntervals();
    updateStartEndTimes();

    // variable declarations and initialization
    double lnProbTimes = 0.0;

    size_t num_extant_sampled = 0;
    size_t num_extant_unsampled = 0;

    double maxb = 0;
    double maxl = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < fbd_taxa.size(); ++i)
    {
        double b = b_i[i];
        double d = d_i[i];

        double o     = fbd_taxa[i].getAgeRange().getMax();
        double y     = fbd_taxa[i].getAgeRange().getMin();
        double o_min = fbd_taxa[i].getMaxAgeRange().getMin();
        double y_max = fbd_taxa[i].getMinAgeRange().getMax();

        size_t bi = l(b);
        size_t di = l(d);
        size_t oi = ages_from_counts ? oldest_intervals[i] : l(o);
        size_t yi = ages_from_counts ? youngest_intervals[i] : l(y);

        // using user-defined first/last age uncertainty
        if ( ages_from_counts == false )
        {
            // check data constraints
            if( oi != l(o_min) || yi != y_max )
            {
                throw("First/last occurrence uncertainty spans multiple intervals");
            }

            if( o_min < y_max )
            {
                throw("First occurrence uncertainty overlaps last occurrence");
            }

            // check model constraints
            if ( !( b > o && ((y == 0.0 && d == 0.0) || (y > 0 && y > d)) && d >= 0.0 ) )
            {
                return RbConstants::Double::neginf;
            }
        }
        // first/last age uncertainty comes from count data
        else
        {
            // check model constraints
            if ( !( b > d && ((y == 0.0 && d == 0.0) || (y > 0 && yi <= di)) && d >= 0.0 ) )
            {
                return RbConstants::Double::neginf;
            }
        }

        // count the number of rho-sampled tips
        num_extant_sampled   += (d == 0.0 && y == 0.0);  // l
        num_extant_unsampled += (d == 0.0 && y > 0.0); // n - m - l

        /*
         * compute speciation/extinction densities
         */

        // find the origin time
        if (b > maxb)
        {
            maxl = birth[bi];
        }

        // include speciation density
        lnProbTimes += log( birth[bi] );

        // multiply by q at the birth time
        lnProbTimes += q(bi, b);

        // include intermediate q terms
        for (size_t j = bi; j < oi; j++)
        {
            lnProbTimes += q_i[j];
        }

        // if there's no uncertainty in o, include factor for the first appearance
        if( ages_from_counts == false && o == o_min )
        {
            lnProbTimes += q(oi, o, true) - q(oi, o);
        }

        // include intermediate q_tilde terms
        for (size_t j = oi; j < di; j++)
        {
            lnProbTimes += q_tilde_i[j];
        }

        // divide by q_tilde at the death time
        lnProbTimes -= q( di, d, true);

        // include extinction density
        if (d > 0.0) lnProbTimes += log( death[di] );

        /*
         * compute sampling density in the first interval
         */
        NaturalNumbersState k_oi = getFossilCount(i,oi);

        // treating missing as positive, i.e. k_oi is always positive
        // (if k_oi is zero, then o == 0.0)
        if ( k_oi.isMissingState() || k_oi.isPositiveState() )
        {
            // include first sampling density
            lnProbTimes += log(fossil[oi]);

            // integrate o over full range of oi
            if ( ages_from_counts == true )
            {
                double delta = std::max(d, times[oi]);
                double delta_plus_Ls_alpha = oi > 0 ? std::min(b, times[oi-1]) : b;

                H[i] = integrateQ(oi, delta_plus_Ls_alpha) - integrateQ(oi, delta);
                lnProbTimes += log( H[i] ) - fossil[oi]*( delta - times[oi] );
            }
            // user-defined uncertainty for o
            else
            {
                if ( o_min == 0.0 )
                {
                    throw("Fossil count > 0 but first occurrence = 0");
                }

                // non-singleton
                if ( oi != yi )
                {
                    // no uncertainty in o
                    if( o == o_min )
                    {
                        lnProbTimes += fossil[oi] * ( o - times[oi] );
                    }
                    // integrate from o_min to o
                    else
                    {
                        // delta = times[oi]
                        H[i] = integrateQ(oi, o) - integrateQ(oi, o_min);
                        lnProbTimes += log( H[i] ); // fossil[oi]*( delta - times[oi] ) = 0
                    }
                }
                // singleton
                else
                {
                    // include last sampling density for extinct species
                    if ( y > 0.0 )
                    {
                        lnProbTimes += log(fossil[oi]);
                    }

                    // no uncertainty in o
                    if( o == o_min )
                    {
                        // no uncertainty in y
                        if ( y == y_max )
                        {
                            lnProbTimes += fossil[oi] * ( o - y );
                        }
                        // integrate over y to y_max
                        else
                        {
                            lnProbTimes += log( exp(-fossil[oi] * y) - exp(-fossil[oi] * y_max) ) - fossil[oi] * o - log(fossil[oi]);
                        }
                    }
                    // integrate from o_min to o
                    else
                    {
                        // delta = y

                        // no uncertainty in y
                        if ( y == y_max )
                        {
                            H[i] = integrateQ(oi, o) - integrateQ(oi, o_min);
                            lnProbTimes += log( H[i] ) - fossil[oi]*( y - times[oi] );
                        }
                        // integrate over y to y_max
                        else
                        {
                            H[i] = integrateQ(oi, o) - integrateQ(oi, o_min);
                            lnProbTimes += log( H[i] );
                            lnProbTimes += log( exp(-fossil[oi] * y) - exp(-fossil[oi] * y_max) ) - fossil[oi] * times[oi] - log(fossil[oi]);
                        }
                    }
                }
            }
        }
        // k_oi fixed
        else
        {
            lnProbTimes += k_oi.getStateIndex() * log(fossil[oi]);

            // integrate o over full range of oi
            if ( ages_from_counts == true )
            {
                double delta = std::max(d, times[oi]);
                double delta_plus_Ls_alpha = oi > 0 ? std::min(b, times[oi-1]) : b;

                H[i] = integrateQ(oi, delta_plus_Ls_alpha, false) - integrateQ(oi, delta, false);
                lnProbTimes += log( H[i] );
            }
            // user-defined uncertainty for o
            else
            {
                if( o_min > 0.0  && k_oi.getStateIndex() == 0 )
                {
                    throw("Fossil count = 0 but first occurrence > 0");
                }
                else if( o == 0.0  && k_oi.getStateIndex() > 0 )
                {
                    throw("Fossil count > 0 but first occurrence = 0");
                }

                // non-singleton
                if ( oi != yi )
                {
                    // integrate from o_min to o
                    if( o != o_min )
                    {
                        H[i] = integrateQ(oi, o, false) - integrateQ(oi, o_min, false);
                        lnProbTimes += log( H[i] );
                    }
                    // no uncertainty in o
                    else
                    {
                        // lnProbTimes += 0;
                    }
                }
                // singleton
                else
                {
                    // integrate from o_min to o
                    if( o != o_min )
                    {
                        // integrate from y to y_max
                        if ( y != y_max )
                        {
                            // ??
                        }
                        // no uncertainty in y
                        else
                        {
                            H[i] = integrateQ(oi, o, false) - integrateQ(oi, o_min, false);
                            lnProbTimes += log( H[i] );
                        }
                    }
                    // integrate from y to y_max
                    else if ( y != y_max )
                    {
                        // ??
                    }
                }
            }
        }

        /*
         * compute sampling density in intervening intervals
         */
        for(size_t j = oi + 1; j < yi; j++)
        {
            NaturalNumbersState k = getFossilCount(i,j);

            double Ls = times[j-1] - times[j];

            // k >= 0
            if ( k.isMissingState() )
            {
                lnProbTimes += fossil[j] * Ls;
            }
            // k > 0
            else if ( k.isPositiveState() )
            {
                lnProbTimes += fossil[j] * Ls + log( 1.0 - exp( - Ls * fossil[j] ) );
            }
            else
            {
                lnProbTimes += k.getStateIndex() * log(fossil[j] * Ls) - RbMath::lnFactorial(k.getStateIndex());
            }
        }

        /*
         * compute sampling density in last interval
         */
        if ( oi != yi )
        {
            NaturalNumbersState k_yi = getFossilCount(i,yi);

            // treat missing as positive, i.e. k_yi is always positive
            // (if k_yi = 0 then o = y = 0.0)
            if ( k_yi.isMissingState() || k_yi.isPositiveState() )
            {
                // include last sampling density for extinct species
                lnProbTimes += log(fossil[yi]);

                // integrate y over full range of yi
                if ( ages_from_counts == true )
                {
                    double Ls = times[yi-1] - std::max(d, times[yi]);

                    lnProbTimes += fossil[yi] * Ls + log( 1.0 - exp( - Ls * fossil[yi] ) );
                }
                // user-defined uncertainty in y
                else
                {
                    // no uncertainty in y
                    if ( y == y_max )
                    {
                        lnProbTimes += fossil[yi] * ( times[yi-1] - y );
                    }
                    // integrate over y to y_max
                    else
                    {
                        lnProbTimes += log( exp(-fossil[oi] * y) - exp(-fossil[oi] * y_max) ) - fossil[oi] * times[yi-1] - log(fossil[oi]);
                    }
                }
            }
            // k_yi fixed
            else
            {
                if( y > 0.0  && k_yi.getStateIndex() == 0 )
                {
                    throw("Fossil count = 0 but last occurrence > 0");
                }

                lnProbTimes += k_yi.getStateIndex() * log(fossil[yi]);

                // integrate y over full range of yi
                if ( ages_from_counts == true )
                {
                    // ??
                }
                // user-defined uncertainty for y
                else if ( y != y_max )
                {
                    // ??
                }
            }
        }

        // abort if likelihood is infinite
        if ( RbMath::isFinite(lnProbTimes) == false )
        {
            return RbConstants::Double::neginf;
        }
    }

    // the origin is not a speciation event
    lnProbTimes -= log(maxl);

    // add the sampled extant tip age term
    if ( homogeneous_rho->getValue() > 0.0)
        lnProbTimes += num_extant_sampled * log( homogeneous_rho->getValue() );
    // add the unsampled extant tip age term
    if ( homogeneous_rho->getValue() < 1.0)
        lnProbTimes += num_extant_unsampled * log( 1.0 - homogeneous_rho->getValue() );

    if ( RbMath::isFinite(lnProbTimes) == false )
    {
        return RbConstants::Double::neginf;
    }

    return lnProbTimes;
}


double AbstractPiecewiseConstantFossilizedRangeProcess::getExtinctionRate( size_t index ) const
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


NaturalNumbersState AbstractPiecewiseConstantFossilizedRangeProcess::getFossilCount( size_t species, size_t interval ) const
{

    // remove the old parameter first
    if ( fossil_count != NULL )
    {
        return NaturalNumbersState(fossil_count->getValue());
    }
    else if( fossil_count_data != NULL )
    {
        interval = ascending ? fossil_count_data->getValue().getNumberOfTaxa() - 1 - interval : interval;

        return fossil_count_data->getValue().getCharacter(species, interval);
    }

    NaturalNumbersState s;
    s.setMissingState(true);

    return s;
}


double AbstractPiecewiseConstantFossilizedRangeProcess::getFossilizationRate( size_t index ) const
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


double AbstractPiecewiseConstantFossilizedRangeProcess::getIntervalTime( size_t index ) const
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


double AbstractPiecewiseConstantFossilizedRangeProcess::getSpeciationRate( size_t index ) const
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
 * \ln\int exp(psi t) q_tilde(t)/q(t) dt
 */
double AbstractPiecewiseConstantFossilizedRangeProcess::integrateQ( size_t i, double t, bool marginalized) const
{
    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == num_intervals - 1 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];

    double diff = b - d - f;
    double bp   = b*f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*bp);
    double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

    double e = exp(-A*dt);

    double diff2 = b + d + (marginalized ? f : -f);
    double tmp = (1+B)/(A-diff2) - e*(1-B)/(A+diff2);
    double intQ = exp(-(diff2-A)*dt/2) * tmp;

//    double tmp = (1.0+B) + e*(1.0-B);
//    double q = 4.0*e / (tmp*tmp);
//    double q_tilde = sqrt(q*exp(-(b+d+f)*dt));
//    double p = b + d + f - A * ((1.0+B)-e*(1.0-B))/((1.0+B)+e*(1.0-B));
//    p /= 2*b;
//    double intQ = q_tilde/q * exp(f*dt) * ( b*p - f)/(b*d-f*d-b*f);

    return intQ;
}


/**
 * return the index i so that t_{i-1} > t >= t_i
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 is origin
 * t_l = 0.0
 */
size_t AbstractPiecewiseConstantFossilizedRangeProcess::l(double t) const
{
    return times.rend() - std::upper_bound( times.rbegin(), times.rend(), t);
}


/**
 * p_i(t)
 */
double AbstractPiecewiseConstantFossilizedRangeProcess::p( size_t i, double t ) const
{
    if (t == 0.0) return 1.0;

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
double AbstractPiecewiseConstantFossilizedRangeProcess::q( size_t i, double t, bool tilde ) const
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
void AbstractPiecewiseConstantFossilizedRangeProcess::updateIntervals() const
{
    std::vector<bool> youngest(fbd_taxa.size(), true);

    for (int i = (int)num_intervals - 1; i >= 0; i--)
    {
        double b = getSpeciationRate(i);
        double d = getExtinctionRate(i);
        double f = getFossilizationRate(i);
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

        if( ages_from_counts == true )
        {
            for(size_t j = 0; j < fbd_taxa.size(); j++)
            {
                NaturalNumbersState s = getFossilCount(i,j);

                if( s.isMissingState() == false && (s.getStateIndex() > 0 || s.isPositiveState()) )
                {
                    oldest_intervals[j] = i;
                    if( youngest[j] )
                    {
                        youngest_intervals[j] = i;
                        youngest[j] = false;
                    }
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
void AbstractPiecewiseConstantFossilizedRangeProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
    else if (oldP == fossil_count)
    {
        fossil_count = static_cast<const TypedDagNode< long >* >( newP );
    }
    else if (oldP == fossil_count_data)
    {
        fossil_count_data = static_cast<const TypedDagNode<HomologousDiscreteCharacterData<NaturalNumbersState> >* >( newP );
    }
}
