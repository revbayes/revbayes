#include "FossilIntervalCounts.h"

#include <cmath>
#include <utility>

#include "RbException.h"
#include "RbConstants.h"
#include "RbMathLogic.h"
#include "TimeInterval.h"

using namespace RevBayesCore;


/**
 * The interval holding an age.
 *
 * Intervals are half-open at their older edge, so an age sitting exactly on a boundary belongs
 * to the younger interval. An age outside the boundaries belongs to no interval, which the
 * callers read as "not sampled".
 */
size_t RevBayesCore::fossilIntervalOf(double age, const std::vector<double> &boundaries, double top)
{
    size_t num_intervals = boundaries.size();

    if ( age < boundaries[0] || age >= top )
    {
        return num_intervals;
    }

    for (size_t j = 0; j + 1 < num_intervals; j++)
    {
        if ( age < boundaries[j+1] )
        {
            return j;
        }
    }

    return num_intervals - 1;
}


void RevBayesCore::checkFossilIntervals(const std::vector<double> &boundaries)
{
    if ( boundaries.size() < 1 )
    {
        throw RbException("The timeline needs at least the present.");
    }

    for (size_t j = 0; j < boundaries.size(); j++)
    {
        if ( RbMath::isFinite( boundaries[j] ) == false || boundaries[j] < 0.0 )
        {
            throw RbException("Interval boundaries must be finite ages.");
        }
        if ( j > 0 && boundaries[j] <= boundaries[j-1] )
        {
            throw RbException("Interval boundaries must be given youngest first, strictly increasing.");
        }
    }
}


bool RevBayesCore::fossilIntervalsEqualWidth(const std::vector<double> &boundaries, double top)
{
    if ( boundaries.size() < 2 )
    {
        return true;
    }

    double first = boundaries[1] - boundaries[0];

    for (size_t j = 1; j < boundaries.size(); j++)
    {
        double width = ( j + 1 < boundaries.size() ? boundaries[j+1] : top ) - boundaries[j];
        if ( RbMath::isFinite( width ) == false )
        {
            continue;
        }
        if ( fabs( width - first ) > 1E-6 * first )
        {
            return false;
        }
    }

    return true;
}


/**
 * Which intervals each taxon was sampled in, from its reported occurrence bins.
 *
 * A bin lying inside one interval counts there. A bin straddling a boundary is counted by the
 * policy: at the interval holding its midpoint, in every interval it touches, or nowhere.
 */
std::vector<std::vector<bool> > RevBayesCore::fossilIntervalPresence(const std::vector<Taxon> &taxa,
                                                                    const std::vector<double> &boundaries,
                                                                    double top,
                                                                    FossilAgeAmbiguity policy)
{
    checkFossilIntervals( boundaries );

    size_t num_intervals = boundaries.size();
    std::vector<std::vector<bool> > presence( taxa.size(), std::vector<bool>( num_intervals, false ) );

    for (size_t i = 0; i < taxa.size(); i++)
    {
        const std::vector<std::pair<TimeInterval, size_t> > &bins = taxa[i].getOccurrences();

        for (size_t b = 0; b < bins.size(); b++)
        {
            double min_age = bins[b].first.getMin();
            double max_age = bins[b].first.getMax();

            size_t youngest = fossilIntervalOf( min_age, boundaries, top );
            size_t oldest   = fossilIntervalOf( max_age, boundaries, top );

            if ( youngest == oldest )
            {
                if ( youngest < num_intervals )
                {
                    presence[i][youngest] = true;
                }
                continue;
            }

            if ( policy == FOSSIL_AGE_MIDPOINT )
            {
                size_t mid = fossilIntervalOf( (min_age + max_age) / 2.0, boundaries, top );
                if ( mid < num_intervals )
                {
                    presence[i][mid] = true;
                }
            }
            else if ( policy == FOSSIL_AGE_OVERLAP )
            {
                size_t from = ( youngest < num_intervals ? youngest : 0 );
                size_t to   = ( oldest   < num_intervals ? oldest   : num_intervals - 1 );
                for (size_t j = from; j <= to; j++)
                {
                    presence[i][j] = true;
                }
            }
        }
    }

    return presence;
}


/**
 * The category counts, all from the presence matrix.
 *
 * A taxon's first appearance is the oldest interval it was sampled in and its last appearance
 * the youngest, which is how the counting methods read a binned record.
 */
RevBayesCore::FossilIntervalCounts RevBayesCore::countFossilIntervals(const std::vector<std::vector<bool> > &presence)
{
    size_t num_intervals = ( presence.empty() ? 0 : presence[0].size() );

    FossilIntervalCounts c;
    c.singleton         = std::vector<double>( num_intervals, 0.0 );
    c.bottom_last       = std::vector<double>( num_intervals, 0.0 );
    c.first_top         = std::vector<double>( num_intervals, 0.0 );
    c.through           = std::vector<double>( num_intervals, 0.0 );
    c.two_timer_older   = std::vector<double>( num_intervals, 0.0 );
    c.two_timer_younger = std::vector<double>( num_intervals, 0.0 );
    c.three_timer       = std::vector<double>( num_intervals, 0.0 );
    c.part_timer        = std::vector<double>( num_intervals, 0.0 );

    for (size_t i = 0; i < presence.size(); i++)
    {
        bool sampled = false;
        size_t first = 0;   // oldest interval sampled in
        size_t last  = 0;   // youngest interval sampled in

        for (size_t j = 0; j < num_intervals; j++)
        {
            if ( presence[i][j] == true )
            {
                if ( sampled == false )
                {
                    last = j;
                    sampled = true;
                }
                first = j;
            }
        }

        if ( sampled == false )
        {
            continue;
        }

        for (size_t j = 0; j < num_intervals; j++)
        {
            bool crosses_older   = ( first > j );
            bool crosses_younger = ( last < j );

            if ( first == j && last == j )
            {
                c.singleton[j]++;
            }
            else if ( crosses_older && last == j )
            {
                c.bottom_last[j]++;
            }
            else if ( first == j && crosses_younger )
            {
                c.first_top[j]++;
            }
            else if ( crosses_older && crosses_younger )
            {
                c.through[j]++;
            }

            // the timer counts need both neighbours, so they skip the two end intervals
            if ( j == 0 || j + 1 == num_intervals )
            {
                continue;
            }

            bool here    = presence[i][j];
            bool older   = presence[i][j+1];
            bool younger = presence[i][j-1];

            if ( here && older )                 c.two_timer_older[j]++;
            if ( here && younger )               c.two_timer_younger[j]++;
            if ( here && older && younger )      c.three_timer[j]++;
            if ( here == false && older && younger ) c.part_timer[j]++;
        }
    }

    return c;
}
