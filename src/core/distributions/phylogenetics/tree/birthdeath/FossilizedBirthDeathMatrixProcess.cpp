#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "DistributionExponential.h"
#include "FossilizedBirthDeathMatrixProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "AbstractFossilizedBirthDeathProcess.h"
#include "MatrixReal.h"
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
 * \param[in]    cdt            Condition of the process (none/survival/#Taxa).
 * \param[in]    tn             Taxa.
 */
FossilizedBirthDeathMatrixProcess::FossilizedBirthDeathMatrixProcess(const DagNode *inspeciation,
                                                                     const DagNode *inextinction,
                                                                     const DagNode *inpsi,
                                                                     const TypedDagNode<double> *inrho,
                                                                     const TypedDagNode< RbVector<double> > *intimes,
                                                                     const std::string &incondition,
                                                                     const std::vector<Taxon> &intaxa,
                                                                     bool complete,
                                                                     bool augmented) :
    TypedDistribution<MatrixReal>(new MatrixReal(intaxa.size(), 2)),
    AbstractFossilizedBirthDeathProcess(inspeciation, inextinction, inpsi, inrho, intimes, intaxa, complete, augmented),
    condition(incondition)
{
    dirty_gamma = std::vector<bool>(fbd_taxa.size(), true);
    gamma_i     = std::vector<size_t>(fbd_taxa.size(), 0);
    gamma_links = std::vector<std::vector<bool> >(fbd_taxa.size(), std::vector<bool>(fbd_taxa.size(), false));

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
FossilizedBirthDeathMatrixProcess* FossilizedBirthDeathMatrixProcess::clone( void ) const
{
    return new FossilizedBirthDeathMatrixProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathMatrixProcess::computeLnProbability( void )
{
    // prepare the probability computation
    updateGamma();

    double lnProb = computeLnProbabilityRanges();

    for( size_t i = 0; i < fbd_taxa.size(); i++ )
    {
        // multiply by the number of possible birth locations
        lnProb += log( gamma_i[i] == 0 ? 1 : gamma_i[i] );
    }

    // condition on survival
    if ( condition == "survival" )
    {
        lnProb -= log( pSurvival(origin,0) );
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
void FossilizedBirthDeathMatrixProcess::updateGamma(bool force)
{
    for (size_t i = 0; i < fbd_taxa.size(); i++)
    {
        if ( dirty_gamma[i] || force )
        {
            double bi = (*this->value)[i][0];
            double di = (*this->value)[i][1];

            if ( force == true ) gamma_i[i] = 0;

            for (size_t j = 0; j < fbd_taxa.size(); j++)
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
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Probability of survival.
 */
double FossilizedBirthDeathMatrixProcess::pSurvival(double start, double end) const
{
    double t = start;

    //std::vector<double> fossil_bak = fossil;

    //std::fill(fossil.begin(), fossil.end(), 0.0);

    double p0 = p(l(t), t);

    //fossil = fossil_bak;

    return 1.0 - p0;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
void FossilizedBirthDeathMatrixProcess::updateStartEndTimes( void )
{
    origin = 0;

    for (size_t i = 0; i < fbd_taxa.size(); i++)
    {
        b_i[i] = (*this->value)[i][0];
        d_i[i] = (*this->value)[i][1];

        origin = std::max(origin, b_i[i]);
    }
}


/**
 * Simulate new speciation times.
 */
void FossilizedBirthDeathMatrixProcess::redrawValue(void)
{
    // incorrect placeholder
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    double max = 0;
    // get the max first occurence
    for (size_t i = 0; i < fbd_taxa.size(); i++)
    {
        double o = fbd_taxa[i].getMaxAge();
        if ( o > max ) max = o;
    }
    
    max *= 1.1;
    
    if (max == 0.0)
    {
        max = 1.0;
    }

    // get random uniform draws
    for (size_t i = 0; i < fbd_taxa.size(); i++)
    {
        double b = fbd_taxa[i].getMaxAge() + rng->uniform01()*(max - fbd_taxa[i].getMaxAge());
        double d = fbd_taxa[i].isExtinct() ? rng->uniform01()*fbd_taxa[i].getMinAge() : 0.0;

        (*this->value)[i][0] = b;
        (*this->value)[i][1] = d;

        redrawOldestOccurrence(i);
    }
}


void FossilizedBirthDeathMatrixProcess::keepSpecialization(DagNode *toucher)
{
    dirty_gamma = std::vector<bool>(fbd_taxa.size(), false);

    AbstractFossilizedBirthDeathProcess::keepSpecialization(toucher);
}


void FossilizedBirthDeathMatrixProcess::restoreSpecialization(DagNode *toucher)
{
    if ( toucher == dag_node )
    {
        updateStartEndTimes();
    }

    AbstractFossilizedBirthDeathProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathMatrixProcess::touchSpecialization(DagNode *toucher, bool touchAll)
{
    if ( toucher == dag_node )
    {
        if ( touched == false )
        {
            stored_likelihood = partial_likelihood;
            stored_o_i = o_i;
            stored_Psi_i = Psi_i;

            std::set<size_t> touched_indices = dag_node->getTouchedElementIndices();

            for ( std::set<size_t>::iterator it = touched_indices.begin(); it != touched_indices.end(); it++)
            {
                size_t i = (*it) / fbd_taxa.size();

                dirty_gamma[i] = true;
                dirty_taxa[i]  = true;

                redrawOldestOccurrence(i);
            }

            updateStartEndTimes();
        }

        touched = true;
    }
    else
    {
        AbstractFossilizedBirthDeathProcess::touchSpecialization(toucher, touchAll);
    }
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FossilizedBirthDeathMatrixProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    AbstractFossilizedBirthDeathProcess::swapParameterInternal(oldP, newP);
}
