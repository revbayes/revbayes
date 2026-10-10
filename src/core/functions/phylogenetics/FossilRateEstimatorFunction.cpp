#include "FossilRateEstimatorFunction.h"

#include <cmath>

#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RlUserInterface.h"

using namespace RevBayesCore;


FossilRateEstimatorFunction::FossilRateEstimatorFunction(const TypedDagNode< RbVector<Taxon> > *t,
                                                         const TypedDagNode< RbVector<double> > *b,
                                                         const TypedDagNode< double > *p,
                                                         const TypedDagNode< double > *mx,
                                                         Method m, FossilAgeAmbiguity a) :
    TypedFunction< MatrixReal >( new MatrixReal( 2, b->getValue().size() + 1 ) ),
    taxa( t ),
    timeline( b ),
    present( p ),
    max_age( mx ),
    method( m ),
    ambiguous( a ),
    warned_max_age( false ),
    warned_widths( false ),
    warned_outside( false )
{
    addParameter( taxa );
    addParameter( timeline );
    addParameter( present );
    addParameter( max_age );

    update();
}


FossilRateEstimatorFunction::~FossilRateEstimatorFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
}


FossilRateEstimatorFunction* FossilRateEstimatorFunction::clone( void ) const
{
    return new FossilRateEstimatorFunction( *this );
}


/**
 * One interval's rate.
 *
 * Foote (2000) eqs 12-13 and 22-23, and Alroy (2008) as given by Warnock et al. (2020) eqs 7-10.
 * An empty denominator, or a zero inside a logarithm, leaves the interval undefined.
 */
double FossilRateEstimatorFunction::rate(const FossilIntervalCounts &c, size_t j, double width, double p_s, bool speciation) const
{
    double nan = RbConstants::Double::nan;

    if ( RbMath::isFinite( width ) == false || width <= 0.0 )
    {
        return nan;
    }

    if ( method == PER_TAXON )
    {
        double total = c.singleton[j] + c.bottom_last[j] + c.first_top[j] + c.through[j];
        if ( total == 0.0 ) return nan;

        double numerator = c.singleton[j] + ( speciation ? c.first_top[j] : c.bottom_last[j] );
        return numerator / total / width;
    }

    if ( method == BOUNDARY_CROSSER )
    {
        double crossing = c.through[j];
        double entering = crossing + ( speciation ? c.first_top[j] : c.bottom_last[j] );
        if ( crossing == 0.0 || entering == 0.0 ) return nan;

        // log(entering/crossing), so an interval where nothing enters reads 0 rather than -0
        return log( entering / crossing ) / width;
    }

    // three-timer
    double three = c.three_timer[j];
    double two = ( speciation ? c.two_timer_younger[j] : c.two_timer_older[j] );
    if ( three == 0.0 || two == 0.0 || p_s <= 0.0 ) return nan;

    return ( log( two / three ) + log( p_s ) ) / width;
}


void FossilRateEstimatorFunction::update( void )
{
    // the process timeline convention: the present, then each rate shift time, youngest first,
    // and the oldest interval unbounded
    const std::vector<double> &shifts = timeline->getValue();
    std::vector<double> boundaries( 1, present->getValue() );
    for (size_t j = 0; j < shifts.size(); j++)
    {
        boundaries.push_back( shifts[j] );
    }

    checkFossilIntervals( boundaries );

    size_t num_intervals = boundaries.size();

    double top = ( max_age != NULL ? max_age->getValue() : RbConstants::Double::inf );

    // A record that stops short of the oldest rate shift leaves that interval no width. Every
    // other undefined case returns NaN, and a column of NaN is easier to act on mid-analysis
    // than an abort, so warn once and let the width test below catch it.
    if ( warned_max_age == false && top <= boundaries[num_intervals-1] )
    {
        RBOUT("Warning: max_age is not older than the oldest rate shift, so that interval has no rate.\n");
        warned_max_age = true;
    }

    // anything older than the oldest edge falls outside every interval and is not counted
    if ( warned_outside == false )
    {
        const RbVector<Taxon> &t = taxa->getValue();
        for (size_t i = 0; i < t.size(); i++)
        {
            if ( t[i].getMaxAge() > top )
            {
                RBOUT("Warning: some occurrences are older than max_age and are not counted.\n");
                warned_outside = true;
                break;
            }
        }
    }

    FossilIntervalCounts c = countFossilIntervals( fossilIntervalPresence( taxa->getValue(), boundaries, top, ambiguous ) );

    // the three-timer sampling probability is pooled over the whole record
    double p_s = RbConstants::Double::nan;
    if ( method == THREE_TIMER )
    {
        double three = 0.0;
        double part = 0.0;
        for (size_t j = 0; j < num_intervals; j++)
        {
            three += c.three_timer[j];
            part += c.part_timer[j];
        }
        p_s = ( three + part > 0.0 ? three / ( three + part ) : RbConstants::Double::nan );

        if ( warned_widths == false && fossilIntervalsEqualWidth( boundaries, top ) == false )
        {
            RBOUT("Warning: the three-timer method assumes intervals of equal length.\n");
            warned_widths = true;
        }
    }

    MatrixReal &v = *value;
    if ( v.getNumberOfRows() != 2 || v.getNumberOfColumns() != num_intervals )
    {
        v = MatrixReal( 2, num_intervals );
    }

    for (size_t j = 0; j < num_intervals; j++)
    {
        double width = ( j + 1 < num_intervals ? boundaries[j+1] : top ) - boundaries[j];

        v[0][j] = rate( c, j, width, p_s, true );
        v[1][j] = rate( c, j, width, p_s, false );
    }
}


void FossilRateEstimatorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if ( oldP == taxa )
    {
        taxa = static_cast<const TypedDagNode< RbVector<Taxon> >* >( newP );
    }
    else if ( oldP == timeline )
    {
        timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if ( oldP == present )
    {
        present = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if ( oldP == max_age )
    {
        max_age = static_cast<const TypedDagNode< double >* >( newP );
    }
}
