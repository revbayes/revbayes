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
#include "RbMathFunctions.h"
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
                                                                                                 bool afc) :
    ascending(false), homogeneous_rho(inrho), timeline( intimes ), fbd_taxa(intaxa), auto_uncertainty(afc), origin(0.0)
{
    dirty_taxa = std::vector<bool>(fbd_taxa.size(), true);
    partial_likelihood = std::vector<double>(fbd_taxa.size(), 0.0);
    stored_likelihood = std::vector<double>(fbd_taxa.size(), 0.0);

    // initialize all the pointers to NULL
    homogeneous_lambda             = NULL;
    homogeneous_mu                 = NULL;
    homogeneous_psi                = NULL;
    heterogeneous_lambda           = NULL;
    heterogeneous_mu               = NULL;
    heterogeneous_psi              = NULL;
    fossil_count_data              = NULL;

    RbException no_timeline_err = RbException("No time intervals provided for piecewise constant fossilized birth death process");

    range_parameters.push_back( homogeneous_rho );

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

    fossil_count_data = dynamic_cast<const TypedDagNode<AbstractHomologousDiscreteCharacterData>* >(incounts);
    range_parameters.push_back( fossil_count_data );

    if ( fossil_count_data != NULL )
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

    b_i = std::vector<double>(fbd_taxa.size(), 0.0);
    d_i = std::vector<double>(fbd_taxa.size(), 0.0);

    lnQ = std::vector<double>(fbd_taxa.size(), 0.0);

    p_i         = std::vector<double>(num_intervals, 1.0);
    q_i         = std::vector<double>(num_intervals, 0.0);
    q_tilde_i   = std::vector<double>(num_intervals, 0.0);

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    updateIntervals();

    // find oldest and youngest intervals
    oldest_intervals = std::vector<size_t>( fbd_taxa.size(), num_intervals - 1 );
    youngest_intervals = std::vector<size_t>( fbd_taxa.size(), num_intervals - 1 );

    for(size_t i = 0; i < fbd_taxa.size(); i++)
    {
        double o     = fbd_taxa[i].getMaxAge();
        double y     = fbd_taxa[i].getMinAge();
        double o_min = fbd_taxa[i].getMaxAgeRange().getMin();
        double y_max = fbd_taxa[i].getMinAgeRange().getMax();

        bool youngest = true;

        for (size_t interval = num_intervals; interval > 0; interval--)
        {
            size_t j = interval - 1;

            if( fossil_count_data != NULL )
            {
                NaturalNumbersState s = getFossilCount(i,j);

                if( s.isMissingState() == false && (s.getStateIndex() > 0 || s.isPositiveState()) )
                {
                    oldest_intervals[i] = j;
                    if( youngest == true )
                    {
                        youngest_intervals[i] = j;
                        youngest = false;
                    }
                }
            }
            else
            {
                if( o >= times[j] )
                {
                    oldest_intervals[i] = j;
                }

                if( y >= times[j] )
                {
                    youngest_intervals[i] = j;
                }
            }
        }

        // check data constraints

        // check consistency between fossil counts and age ranges
        if( fossil_count_data != NULL && (oldest_intervals[i] != l(o) || youngest_intervals[i] != l(y)) )
        {
            std::stringstream ss;
            ss << "Inconsistent age range and fossil count data for taxon " << fbd_taxa[i].getName();
            throw(RbException(ss.str()));
        }

        // check data constraints
        if( oldest_intervals[i] != l(o_min) || youngest_intervals[i] != l(y_max) )
        {
            std::stringstream ss;
            ss << "First/last occurrence uncertainty spans multiple intervals for taxon " << fbd_taxa[i].getName();
            throw(RbException(ss.str()));
        }

        if( o_min < y_max )
        {
            std::stringstream ss;
            ss << "First occurrence uncertainty overlaps last occurrence for taxon " << fbd_taxa[i].getName();
            throw(RbException(ss.str()));
        }
    }
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

    // add the fossil tip age terms
    for (size_t i = 0; i < fbd_taxa.size(); ++i)
    {
        double b = b_i[i];
        double d = d_i[i];

        double o     = fbd_taxa[i].getMaxAge();
        double y     = fbd_taxa[i].getMinAge();
        double o_min = fbd_taxa[i].getMaxAgeRange().getMin();
        double y_max = fbd_taxa[i].getMinAgeRange().getMax();

        // count the number of rho-sampled tips
        num_extant_sampled   += (d == 0.0 && y == 0.0);  // l
        num_extant_unsampled += (d == 0.0 && y > 0.0); // n - m - l

        if ( dirty_taxa[i] == true )
        {
            partial_likelihood[i] = 0.0;

            size_t bi = l(b);
            size_t di = l(d);
            size_t oi = oldest_intervals[i];
            size_t yi = youngest_intervals[i];

            // check model constraints
            if ( !( b > o && ((y == 0.0 && d == 0.0) || (y > 0 && y > d)) && d >= 0.0 ) )
            {
                return RbConstants::Double::neginf;
            }

            // include speciation density
            partial_likelihood[i] += log( birth[bi] );

            /*
             * compute speciation/extinction densities
             */

            // multiply by q at the birth time
            partial_likelihood[i] += q(bi, b);

            // include intermediate q terms
            for (size_t j = bi; j < oi; j++)
            {
                partial_likelihood[i] += q_i[j];
            }

            // if there's no uncertainty in o, include factor for the first appearance
            if( auto_uncertainty == false && o == o_min )
            {
                partial_likelihood[i] += q(oi, o, true) - q(oi, o);
            }

            // include intermediate q_tilde terms
            for (size_t j = oi; j < di; j++)
            {
                partial_likelihood[i] += q_tilde_i[j];
            }

            // divide by q_tilde at the death time
            partial_likelihood[i] -= q( di, d, true);

            // include extinction density
            if (d > 0.0) partial_likelihood[i] += log( death[di] );

            /*
             * compute sampling density in the first interval
             */
            NaturalNumbersState k_oi = getFossilCount(i,oi);

            // make sure extant taxa don't have fossils
            if ( o == 0.0 )
            {
                if ( fossil_count_data != NULL && (k_oi.isAmbiguous() == true || k_oi.getStateIndex() > 0) )
                {
                    std::stringstream ss;
                    ss << "First occurrence fossil count > 0 for extant taxon " << fbd_taxa[i].getName();
                    throw(RbException(ss.str()));
                }
            }
            // treating missing as positive, i.e. k_oi is always positive
            // (if k_oi is zero, then o == 0.0)
            else if ( k_oi.isMissingState() || k_oi.isPositiveState() )
            {
                // include first sampling density
                partial_likelihood[i] += log(fossil[oi]);

                // integrate o over full range of oi
                if ( auto_uncertainty == true )
                {
                    // non-singleton
                    // (Case 2)
                    if ( oi != yi )
                    {
                        double delta = times[oi];
                        double delta_plus_Ls = oi > 0 ? std::min(b, times[oi-1]) : b;

                        lnQ[i] = H(oi,times[oi],delta_plus_Ls) - H(oi,times[oi],delta);
                        lnQ[i] -= Z(0,oi,times[oi],delta_plus_Ls) - Z(0,oi,times[oi],delta);
                        lnQ[i] = log(lnQ[i]);
                    }
                    // singleton
                    // (Case 6)
                    else
                    {
                        double delta = std::max(d, times[oi]);
                        double delta_plus_Ls = oi > 0 ? std::min(b, times[oi-1]) : b;

                        lnQ[i] = H(oi,delta,delta_plus_Ls) - H(oi,delta_plus_Ls,delta_plus_Ls) - (H(oi,delta,delta) - H(oi,delta_plus_Ls,delta));
                        lnQ[i] /= fossil[oi];
                        lnQ[i] -= Z(1,oi,delta,delta_plus_Ls) - Z(1,oi,delta,delta);
                        lnQ[i] = log(lnQ[i]);
                    }

                    partial_likelihood[i] += lnQ[i];
                }
                // user-defined uncertainty for o
                else
                {
                    // non-singleton
                    if ( oi != yi )
                    {
                        // no uncertainty in o
                        // (Eq 6)
                        if( o == o_min )
                        {
                            partial_likelihood[i] += fossil[oi] * ( o - times[oi] );
                        }
                        // integrate from o_min to o
                        // (Case 2)
                        else
                        {
                            lnQ[i] = H(oi,times[oi],o) - H(oi,times[oi],o_min);
                            lnQ[i] -= Z(0,oi,times[oi],o) - Z(0,oi,times[oi],o_min);
                            lnQ[i] = log(lnQ[i]);

                            partial_likelihood[i] += lnQ[i];
                        }
                    }
                    // singleton
                    else
                    {
                        // include last sampling density for extinct species
                        if ( o_min != y_max && y > 0.0 )
                        {
                            partial_likelihood[i] += log(fossil[oi]);
                        }

                        // no uncertainty in o
                        if ( o == o_min )
                        {
                            // no uncertainty in y
                            // (Eq 6)
                            if ( y == y_max )
                            {
                                partial_likelihood[i] += fossil[oi] * ( o - y );
                            }
                            // integrate over y to y_max
                            // (Case 6)
                            else
                            {
                                lnQ[i] = H(oi,y_max,o) - H(oi,y,o);
                                lnQ[i] /= fossil[oi];
                                lnQ[i] += Z(1,oi,y_max,o) - Z(1,oi,y,o) - (Z(1,oi,y_max,o_min) - Z(1,oi,y,o_min));
                                lnQ[i] = log(lnQ[i]);
                                partial_likelihood[i] += lnQ[i];
                            }
                        }
                        // integrate from o_min to o
                        else
                        {
                            // delta = y

                            // no uncertainty in y
                            // (Case 6)
                            if ( y == y_max )
                            {
                                lnQ[i] = H(oi,y,o) - H(oi,y,o_min);
                                lnQ[i] += Z(0,oi,y,o) - Z(0,oi,y,o_min);
                            }
                            // integrate over y to y_max
                            // (Case 6)
                            else
                            {
                                lnQ[i] = H(oi,y_max,o) - H(oi,y_max,o_min) - H(oi,y,o) + H(oi,y,o_min);
                                lnQ[i] /= fossil[oi];
                                lnQ[i] += Z(1,oi,y_max,o) - Z(1,oi,y,o) - (Z(1,oi,y_max,o_min) - Z(1,oi,y,o_min));
                            }

                            lnQ[i] = log(lnQ[i]);
                            partial_likelihood[i] += lnQ[i];
                        }
                    }
                }
            }
            // k_oi fixed
            else
            {
                size_t k = k_oi.getStateIndex();

                partial_likelihood[i] += k * log(fossil[oi]);

                // integrate o over full range of oi
                if ( auto_uncertainty == true )
                {
                    // non-singleton
                    // (Case 1)
                    if ( oi != yi )
                    {
                        double delta = std::max(d, times[oi]);
                        double delta_plus_Ls = oi > 0 ? std::min(b, times[oi-1]) : b;

                        lnQ[i] = Z(k, oi, times[oi], delta_plus_Ls) - Z(k, oi, times[oi], delta);
                        lnQ[i] = log(lnQ[i]);
                        partial_likelihood[i] += lnQ[i] - RbMath::lnFactorial(k);
                    }
                    // singleton
                    // (Case 5)
                    else
                    {
                        double delta = std::max(d, times[oi]);
                        double delta_plus_Ls = oi > 0 ? std::min(b, times[oi-1]) : b;

                        lnQ[i] = Z(k+1,oi,delta,delta_plus_Ls) - Z(k+1,oi,delta,delta);
                        lnQ[i] = log(lnQ[i]);
                        partial_likelihood[i] += lnQ[i] - RbMath::lnFactorial(k+1);
                    }
                }
                // user-defined uncertainty for o
                else
                {
                    // non-singleton
                    if ( oi != yi )
                    {
                        // integrate from o_min to o
                        // (Case 1)
                        if ( o != o_min )
                        {
                            lnQ[i] = Z(k, oi, times[oi], o) - Z(k, oi, times[oi], o_min);
                            lnQ[i] = log(lnQ[i]);
                            partial_likelihood[i] += lnQ[i] - RbMath::lnFactorial(k);
                        }
                        // no uncertainty in o
                        // (Eq 5)
                        else
                        {
                            partial_likelihood[i] += k * log(o - times[oi]) - RbMath::lnFactorial(k);
                        }
                    }
                    // singleton
                    else
                    {
                        // integrate from o_min to o
                        if ( o != o_min )
                        {
                            // integrate from y to y_max
                            // (Case 5)
                            if ( y != y_max )
                            {
                                lnQ[i] = Z(k+1, oi, y_max, o) - Z(k+1, oi, y, o) - (Z(k+1, oi, y_max, o_min) - Z(k+1, oi, y, o_min));
                                lnQ[i] = log(-lnQ[i]);
                                partial_likelihood[i] += lnQ[i] - RbMath::lnFactorial(k+1);
                            }
                            // no uncertainty in y
                            // (Case 1)
                            else
                            {
                                lnQ[i] = Z(k, oi, y, o) - Z(k, oi, y, o_min);
                                lnQ[i] = log(lnQ[i]);
                                partial_likelihood[i] += lnQ[i] - RbMath::lnFactorial(k);
                            }
                        }
                        // no uncertainty in o
                        else
                        {
                            // integrate from y to y_max
                            // (Case 3)
                            if ( y != y_max )
                            {
                                partial_likelihood[i] += log(pow(o-y,k+1) - pow(o-y_max,k+1)) - RbMath::lnFactorial(k + 1);
                            }
                            // no uncertainty in y
                            // (Eq 5)
                            else if( o != y )
                            {
                                partial_likelihood[i] += k * log(o - y) - RbMath::lnFactorial(k);
                            }
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
                // (Eq 6)
                if ( k.isMissingState() )
                {
                    partial_likelihood[i] += fossil[j] * Ls;
                }
                // k > 0
                // (Eq 6)
                else if ( k.isPositiveState() )
                {
                    partial_likelihood[i] += fossil[j] * Ls + log( 1.0 - exp( - Ls * fossil[j] ) );
                }
                // k fixed
                // (Eq 5)
                else
                {
                    partial_likelihood[i] += k.getStateIndex() * log(fossil[j] * Ls) - RbMath::lnFactorial(k.getStateIndex());
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
                    if ( y > 0.0 )
                    {
                        partial_likelihood[i] += log(fossil[yi]);
                    }

                    // integrate y over full range of yi
                    // (Case 4)
                    if ( auto_uncertainty == true )
                    {
                        double Ls = times[yi-1] - std::max(d, times[yi]);

                        double tmp = fossil[yi] * (1.0 - exp(fossil[yi] * Ls));
                        tmp += Ls;
                        partial_likelihood[i] += log(tmp);
                    }
                    // user-defined uncertainty in y
                    else
                    {
                        // no uncertainty in y
                        // (Eq 6)
                        if ( y == y_max )
                        {
                            double tmp = exp(fossil[yi] * ( times[yi-1] - y ));
                            tmp = expm1(fossil[yi] * ( times[yi-1] - y ));
                            partial_likelihood[i] += log(tmp);
                        }
                        // integrate over y to y_max
                        // (Case 4)
                        else
                        {
                            double tmp = fossil[yi] *( exp(fossil[yi] * (times[yi-1]-y)) - exp(fossil[yi] * (times[yi-1]-y_max)) );
                            tmp += (times[yi-1]-y_max) - (times[yi-1]-y);
                            partial_likelihood[i] += log(tmp);
                        }
                    }
                }
                // k_yi fixed
                else
                {
                    size_t k = k_yi.getStateIndex();

                    partial_likelihood[i] += k * log(fossil[yi]);

                    // integrate y over full range of yi
                    // (Case 3)
                    if ( auto_uncertainty == true )
                    {
                        double Ls = times[yi-1] - std::max(d, times[yi]);
                        partial_likelihood[i] += (k+1)*log(Ls) - RbMath::lnFactorial(k + 1);
                    }
                    // user-defined uncertainty for y
                    // (Case 3)
                    else if ( y != y_max )
                    {
                        partial_likelihood[i] += log(pow(times[yi-1]-y,k+1) - pow(times[yi-1]-y_max,k+1)) - RbMath::lnFactorial(k + 1);
                    }
                    // no uncertainty in y
                    // (Eq 5)
                    else
                    {
                        partial_likelihood[i] += k * log(times[yi-1]-y_max) - RbMath::lnFactorial(k);
                    }
                }
            }
        }

        lnProbTimes += partial_likelihood[i];
    }

    // (the origin is not a speciation event)
    lnProbTimes -= log( birth[l(origin)] );

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
    if( fossil_count_data != NULL )
    {
        interval = ascending ? num_intervals - 1 - interval : interval;

        const NaturalNumbersState& count = dynamic_cast<const NaturalNumbersState&>(fossil_count_data->getValue().getCharacter(species, interval));
        return count;
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
 * \int exp(psi (t-t_i)) q_tilde(t)/q(t) dt
 */
double AbstractPiecewiseConstantFossilizedRangeProcess::H(size_t i, double x, double t) const
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

    double diff2 = b + d -f;
    double tmp = (1+B)/(A-diff2) - e*(1-B)/(A+diff2);
    double H = exp(-f*(x-ti) ) * exp(-(diff2-A)*dt/2) * tmp;

    return H;
}


/**
 * \int (t-t_j)^k q_tilde(t)/q(t) dt
 */
double AbstractPiecewiseConstantFossilizedRangeProcess::Z(size_t k, size_t i, double x, double t) const
{
    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == num_intervals - 1 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];

    double diff = b - d - f;
    double bp   = b*f;

    double A = sqrt( diff*diff + 4.0*bp);
    double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

    double sum = b + d + f;
    double alpha = 0.5*(A+sum);

    double tmp1 = pow(-2,k) * (1+B) * exp(-(alpha-A)*(x-ti)) / pow(A-sum,k+1);
    double tmp2 = pow(2,k)  * (1-B) * exp(-alpha*(x-ti))     / pow(A+sum,k+1);

    try
    {
        tmp1 *= RbMath::incompleteGamma((alpha-A)*(t-x), k+1, false, false);
        tmp2 *= RbMath::incompleteGamma(alpha*(t-x), k+1, false, false);
    }
    catch(RbException&)
    {
        return RbConstants::Double::nan;
    }

    return tmp1 - tmp2;
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
    for (size_t interval = num_intervals; interval > 0; interval--)
    {
        size_t i = interval - 1;

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
    else if (oldP == fossil_count_data)
    {
        fossil_count_data = static_cast<const TypedDagNode<AbstractHomologousDiscreteCharacterData>* >( newP );
    }
}
