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
                                                                     size_t truncate_at,
                                                                                                                                          const TypedDagNode<double> *inorigin,
                                                                     TypedDistribution<double> *inoriginprior,
                                                                     bool report_int) :
    TypedDistribution<MatrixReal>(new MatrixReal(intaxa.size(), 2)),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete_record, truncate_at, inorigin, inoriginprior)
{
    report_internally = report_int;

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
    if ( tf == true ) warnIfNoResampleMove();
    if ( tf == true ) warnIfNoReportingNode();
}

FossilizedBirthDeathRangeProcess* FossilizedBirthDeathRangeProcess::clone( void ) const
{
    return new FossilizedBirthDeathRangeProcess( *this );
}


/**
 * Set the matrix value (e.g. clamping to fixed birth/death ages). A clamp replaces the
 * b/d that redrawValue drew the augmented ages against, so re-clip any age now out of
 * range -- otherwise a clamped chain can start at lnProb = -inf. MCMC moves edit the
 * value in place instead, leaving out-of-range ages for the constraints to reject.
 */
void FossilizedBirthDeathRangeProcess::setValue(MatrixReal *v, bool force)
{
    TypedDistribution<MatrixReal>::setValue(v, force);

    for (size_t i = 0; i < taxa.size(); i++)
    {
        double d  = (*this->value)[i][1];
        double b  = (*this->value)[i][0];

        // oldest age: valid range [max(o_i,d), min(max_age,b))
        double lo = std::max( o_i[i], d );
        double hi = std::min( taxa[i].getMaxAge(), b );
        if ( hi > lo && ( first[i] < lo || first[i] >= b ) )
        {
            first[i] = 0.5 * ( lo + hi );
        }

        // youngest age: valid range [max(d,min_age), min(first,y_i)]
        double lo_y = std::max( d, taxa[i].getMinAge() );
        double hi_y = std::min( first[i], y_i[i] );
        if ( hi_y > lo_y && ( last[i] < lo_y || last[i] > hi_y ) )
        {
            last[i] = 0.5 * ( lo_y + hi_y );
        }
    }
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
            double bi = (*this->value)[i][0];
            double di = (*this->value)[i][1];

            if ( force == true ) gamma_i[i] = 0;

            for (size_t j = 0; j < taxa.size(); j++)
            {
                if (i == j) continue;

                double bj = (*this->value)[j][0];
                double dj = (*this->value)[j][1];

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
void FossilizedBirthDeathRangeProcess::updateStartEndTimes( void )
{
    max_birth = 0;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        b_i[i] = (*this->value)[i][0];
        d_i[i] = (*this->value)[i][1];

        if ( b_i[i] > b_i[max_birth] ) max_birth = i;
    }

    if ( origin_age != NULL )
    {
        origin = origin_age->getValue();
    }
    else
    {
        origin = b_i[max_birth];
    }
}


/**
 * Simulate new speciation times.
 */
void FossilizedBirthDeathRangeProcess::redrawValue(void)
{
    // draw an initial range per taxon (the tree process shares this and hangs a topology on it);
    // updateStartEndTimes overwrites b_i/d_i once the matrix is set
    drawRanges();

    for (size_t i = 0; i < taxa.size(); i++)
    {
        (*this->value)[i][0] = b_i[i];
        (*this->value)[i][1] = d_i[i];
    }
}


void FossilizedBirthDeathRangeProcess::keepSpecialization(const DagNode *toucher)
{
    dirty_gamma = std::vector<bool>(taxa.size(), false);

    AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(toucher);
}

void FossilizedBirthDeathRangeProcess::restoreSpecialization(const DagNode *toucher)
{
    AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathRangeProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    if ( toucher == dag_node )
    {
        if ( touched == false )
        {
            stored_likelihood = partial_likelihood;
            stored_Psi = Psi;

            std::set<size_t> touched_indices = dag_node->getTouchedElementIndices();

            for ( std::set<size_t>::iterator it = touched_indices.begin(); it != touched_indices.end(); it++)
            {
                size_t i = (*it) / 2; // (birth,death) matrix is N x 2, row-major: linear index / 2 = taxon

                dirty_gamma[i] = true;
                dirty_psi[i]   = true;
                dirty_taxa[i]  = true;

            }

        }

        touched = true;
    }
    else
    {
        AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(toucher, touchAll);
    }
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
