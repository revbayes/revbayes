#include "BirthDeathWithRateshifts.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "RbConstants.h"
#include "RbMathLogic.h"
#include "Taxon.h"

using namespace RevBayesCore;


BirthDeathWithRateshifts::BirthDeathWithRateshifts(const DagNode *inspeciation,
                                                   const DagNode *inextinction,
                                                   const DagNode *inpsi,
                                                   const TypedDagNode<double> *inrho,
                                                   const TypedDagNode< RbVector<double> > *intimes,
                                                   const std::string &incondition,
                                                   const std::vector<Taxon> &intaxa,
                                                   const std::string &reporting,
                                                   size_t truncate_at,
                                                   bool resample,
                                                   const TypedDagNode<double> *inorigin,
                                                   bool report_int) :
    FossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, reporting, truncate_at, resample, inorigin, NULL, report_int)
{

}


BirthDeathWithRateshifts* BirthDeathWithRateshifts::clone( void ) const
{
    return new BirthDeathWithRateshifts( *this );
}


/**
 * BDS (Silvestro et al. 2019) treats lineages as independent: no coexistence (gamma) factor, and
 * each range is normalized by the fossil non-detection probability e^{-psi_b_d} over [d,b].
 */
double BirthDeathWithRateshifts::computeLnProbability( void )
{
    prepareProbComputation();

    double lnProb = 0.0;

    size_t num_rho_sampled = 0;

    // add the fossil tip age terms
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        double b = this->getValue()[i][0];
        double d = this->getValue()[i][1];

        double max_age = taxa[i].getMaxAge();

        double present = times.front();

        // check model constraints
        if ( !( b > o_i[i] && b > d && y_i[i] >= d && d >= present ) )
        {
            return RbConstants::Double::neginf;
        }
        if ( (d > present) != taxa[i].isExtinct() )
        {
            return RbConstants::Double::neginf;
        }

        // count the number of rho-sampled tips (see computeLnProbabilityRanges)
        num_rho_sampled += (d == present);

        if ( dirty_taxa[i] == true )
        {
            size_t bi = findIndex(b);
            size_t di = findIndex(d);

            partial_likelihood[i] = 0.0;

            // include speciation density
            partial_likelihood[i] += log( birth[bi] );

            // skip the rest for extant taxa with no fossil samples
            if ( max_age == present )
            {
                continue;
            }

            // include extinction density
            if (d > present) partial_likelihood[i] += log( death[di] );

            double psi_b_d = 0.0;

            // include poisson density
            for ( size_t j = di; j <= bi; j++ )
            {
                double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

                double dt = std::min(b, t_0) - std::max(d, times[j]);

                partial_likelihood[i] -= (birth[j] + death[j])*dt;

                psi_b_d += fossil[j]*dt;
            }

            // fossil non-sampling normalization over [d,b] (range-process term)
            partial_likelihood[i] -= psi_b_d;

            if ( condition == "sampling" )
            {
                partial_likelihood[i] -= log(-expm1(-psi_b_d));
            }

            // inline reporting term, unless a dnFossilRecord node supplies it
            if ( report_internally && dirty_psi[i] )
            {
                b_i[i] = b;
                d_i[i] = d;
                Psi[i] = computeLnFossilRecord(i);
                if ( Psi[i] == RbConstants::Double::neginf )
                {
                    return RbConstants::Double::neginf;
                }
            }

            if ( report_internally ) partial_likelihood[i] += Psi[i];
        }

        lnProb += partial_likelihood[i];
    }

    // add the sampled extant tip age term
    if ( homogeneous_rho->getValue() > 0.0)
    {
        lnProb += num_rho_sampled * log( homogeneous_rho->getValue() );
    }

    if ( RbMath::isFinite(lnProb) == false )
    {
        return RbConstants::Double::neginf;
    }

    return lnProb;
}
