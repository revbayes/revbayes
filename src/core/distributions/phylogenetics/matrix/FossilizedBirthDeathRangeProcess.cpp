#include "FossilizedBirthDeathRangeProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "DistributionExponential.h"
#include "MatrixReal.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "RbMathFunctions.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "TimeInterval.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

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
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 */
FossilizedBirthDeathRangeProcess::FossilizedBirthDeathRangeProcess(const DagNode *inspeciation,
                                                                     const DagNode *inextinction,
                                                                     const DagNode *inpsi,
                                                                     const TypedDagNode<double> *inrho,
                                                                     const TypedDagNode< RbVector<double> > *intimes,
                                                                     const std::string &incondition,
                                                                     const std::vector<Taxon> &intaxa,
                                                                     bool complete_record,
                                                                                                                                          const TypedDagNode<double> *inorigin,
                                                                     TypedDistribution<double> *inoriginprior,
                                                                     bool insurvivors,
                                                                     double inpresent) :
    TypedDistribution<MatrixReal>(new MatrixReal(intaxa.size(), 4)),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete_record, inorigin, inoriginprior, insurvivors, inpresent)
{

    dirty_gamma = std::vector<bool>(taxa.size(), true);
    gamma_i     = std::vector<size_t>(taxa.size(), 0);
    gamma_links = std::vector<std::vector<bool> >(taxa.size(), std::vector<bool>(taxa.size(), false));

    for(std::vector<const DagNode*>::iterator it = range_parameters.begin(); it != range_parameters.end(); it++)
    {
        addParameter(*it);
    }

    redrawValue();
    updateGamma(true);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */

void FossilizedBirthDeathRangeProcess::setMcmcMode(bool tf)
{
    TypedDistribution<MatrixReal>::setMcmcMode(tf);
    // the appearances live in the value here, where the element moves reach them, so only
    // the tree process needs to warn that nothing samples them
    if ( tf == true ) warnIfNoReportingNode();
}

FossilizedBirthDeathRangeProcess* FossilizedBirthDeathRangeProcess::clone( void ) const
{
    return new FossilizedBirthDeathRangeProcess( *this );
}


/**
 * Set the matrix value (e.g. clamping to fixed birth/death ages). A clamp replaces the
 * b/d that redrawValue drew the appearances against, so re-clip any age now out of
 * range -- otherwise a clamped chain can start at lnProb = -inf. MCMC moves edit the
 * value in place instead, leaving out-of-range ages for the constraints to reject.
 */
void FossilizedBirthDeathRangeProcess::setValue(MatrixReal *v, bool force)
{
    TypedDistribution<MatrixReal>::setValue(v, force);

    for (size_t i = 0; i < taxa.size(); i++)
    {
        double d  = (*this->value)[i][0];
        double b  = (*this->value)[i][3];

        // read them from the value, not the table, which still holds the draw this replaces
        ranges[i].last  = (*this->value)[i][1];
        ranges[i].first = (*this->value)[i][2];

        // a value from outside carries no status, so take the one its death implies
        ranges[i].survived = ( d == getPresent() );

        // oldest age: valid range [max(first_min,d), min(max_age,b))
        double lo = std::max( ranges[i].first_min, d );
        double hi = std::min( taxa[i].getMaxAge(), b );
        if ( hi > lo && ( ranges[i].first < lo || ranges[i].first >= b ) )
        {
            ranges[i].first = 0.5 * ( lo + hi );
        }

        // youngest age: valid range [max(d,min_age), min(first,last_max)]
        double lo_y = std::max( d, taxa[i].getMinAge() );
        double hi_y = std::min( ranges[i].first, ranges[i].last_max );
        if ( hi_y > lo_y && ( ranges[i].last < lo_y || ranges[i].last > hi_y ) )
        {
            ranges[i].last = 0.5 * ( lo_y + hi_y );
        }

        // the value is what the density reads, so a clip that only moved the members is lost
        (*this->value)[i][1] = ranges[i].last;
        (*this->value)[i][2] = ranges[i].first;
    }

    repairRanges();

    updateRanges();
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathRangeProcess::computeLnProbability( void )
{
    updateGamma();

    double lnProb = computeLnProbabilityRanges();

    for( size_t i = 0; i < taxa.size(); i++ )
    {
        if ( gamma_i[i] == 0 )
        {
            // only the origin buds from nothing
            if ( i == max_birth ) continue;

            return RbConstants::Double::neginf;
        }

        // multiply by the number of possible birth locations
        lnProb += log( gamma_i[i] );
    }

    return lnProb;
}


/**
 * Compute the number of ranges that intersect with range i
 *
 * \param[in]    i      index of range for which to compute gamma
 *
 * \return Small gamma
 */
void FossilizedBirthDeathRangeProcess::updateGamma(bool force)
{
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( dirty_gamma[i] || force )
        {
            double bi = (*this->value)[i][3];
            double di = (*this->value)[i][0];

            if ( force == true ) gamma_i[i] = 0;

            for (size_t j = 0; j < taxa.size(); j++)
            {
                if (i == j) continue;

                double bj = (*this->value)[j][3];
                double dj = (*this->value)[j][0];

                bool linki = ( bi < bj && bi > dj );
                bool linkj = ( bj < bi && bj > di );

                if ( gamma_links[i][j] != linki && force == false )
                {
                    gamma_i[i] += linki ? 1 : -1;
                }
                if ( gamma_links[j][i] != linkj && force == false )
                {
                    gamma_i[j] += linkj ? 1 : -1;
                }

                if ( force == true ) gamma_i[i] += linki;

                gamma_links[i][j] = linki;
                gamma_links[j][i] = linkj;
            }
        }
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
void FossilizedBirthDeathRangeProcess::updateRanges( void )
{
    max_birth = 0;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        ranges[i].birth = (*this->value)[i][3];
        ranges[i].death = (*this->value)[i][0];
        ranges[i].first = (*this->value)[i][2];
        ranges[i].last  = (*this->value)[i][1];

        if ( ranges[i].birth > ranges[max_birth].birth ) max_birth = i;
    }

    if ( origin_age != NULL )
    {
        origin = origin_age->getValue();
    }
    else
    {
        origin = ranges[max_birth].birth;
    }
}


/**
 * Simulate new speciation times.
 */
void FossilizedBirthDeathRangeProcess::redrawValue(void)
{
    // draw an initial range per taxon (the tree process shares this and hangs a topology on it);
    // updateRanges overwrites range_start/range_end once the matrix is set
    drawRanges();

    for (size_t i = 0; i < taxa.size(); i++)
    {
        (*this->value)[i][0] = ranges[i].death;
        (*this->value)[i][1] = ranges[i].last;
        (*this->value)[i][2] = ranges[i].first;
        (*this->value)[i][3] = ranges[i].birth;
    }

    repairRanges();

    // repairRanges writes the value, and the value is what the density reads through the table
    updateRanges();
}


/**
 * A taxon reporting fewer than two occurrences has one appearance, not two: tau_1 and tau_K are
 * the same quantity, and the value carries a redundant copy. A move that slides the two columns
 * apart is repaired here rather than rejected, which is forced (the reverse is equally forced, so
 * the ratio is 1) and costs no mixing. Taxa with both extremes are ordered rather than equal, and
 * the density already rejects those out of order.
 *
 * On touch, not on score: the move stored the value, so a rejected proposal takes the repair with
 * it. A density that writes to the value it scores leaves the write behind.
 */
void FossilizedBirthDeathRangeProcess::repairRanges( const std::set<size_t> &touched )
{
    for ( std::set<size_t>::const_iterator it = touched.begin(); it != touched.end(); it++ )
    {
        size_t i = (*it) / 4;
        size_t c = (*it) % 4;

        if ( ranges[i].singleton == false ) continue;

        // follow whichever column the move wrote; the sibling is stored here, since the move
        // stored only the element it touched
        size_t sib;
        if      ( c == 1 ) sib = 2;
        else if ( c == 2 ) sib = 1;
        else continue;

        stored_repairs.push_back( std::make_pair( i*4 + sib, (*this->value)[i][sib] ) );
        (*this->value)[i][sib] = (*this->value)[i][c];
    }
}


void FossilizedBirthDeathRangeProcess::repairRanges( void )
{
    // no move to follow (a fresh draw, or a value set from outside): tau_1 is the reported one
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( ranges[i].singleton ) (*this->value)[i][1] = (*this->value)[i][2];
    }
}


void FossilizedBirthDeathRangeProcess::keepSpecialization(const DagNode *toucher)
{
    stored_repairs.clear();

    dirty_gamma = std::vector<bool>(taxa.size(), false);

    AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(toucher);
}

void FossilizedBirthDeathRangeProcess::restoreSpecialization(const DagNode *toucher)
{
    // undo the repair's writes newest first, so repeated writes to one element unwind correctly
    for ( std::vector<std::pair<size_t,double> >::reverse_iterator it = stored_repairs.rbegin(); it != stored_repairs.rend(); ++it )
    {
        (*this->value)[ it->first / 4 ][ it->first % 4 ] = it->second;
    }
    stored_repairs.clear();

    AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathRangeProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    if ( toucher == dag_node )
    {
        if ( touched == false )
        {
            stored_likelihood = partial_likelihood;
            stored_ranges     = ranges;

            std::set<size_t> touched_indices = dag_node->getTouchedElementIndices();

            for ( std::set<size_t>::iterator it = touched_indices.begin(); it != touched_indices.end(); it++)
            {
                size_t i = (*it) / 4; // N x 4 row-major, so the linear index over the columns is the taxon

                dirty_gamma[i] = true;
                dirty_taxa[i]  = true;

            }

            repairRanges( touched_indices );
        }

        touched = true;
    }
    else
    {
        AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(toucher, touchAll);
    }

    // the proposal has already written the value, so pull now
    updateRanges();
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FossilizedBirthDeathRangeProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(oldP, newP);
}
