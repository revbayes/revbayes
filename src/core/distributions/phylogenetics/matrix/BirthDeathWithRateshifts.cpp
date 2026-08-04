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
                                                   bool complete_record,
                                                   const TypedDagNode<double> *inorigin,
                                                   double inpresent) :
    FossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete_record, inorigin, NULL, true, inpresent)
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
/**
 * The likelihood is the range process's, so only the birth-death term below differs. Nothing
 * about rho does: the status flag is data here too, and a taxon reported extinct is free to
 * have survived unseen and pay 1-rho.
 */
double BirthDeathWithRateshifts::computeLnProbability( void )
{
    // no gamma factor: lineages are independent under complete lineage sampling
    return computeLnProbabilityRanges();
}


/**
 * Birth-death density of range i under Silvestro et al. (2019): exponential waiting over the
 * range and a Poisson record on it, with no survival probability p(t) anywhere.
 */
double BirthDeathWithRateshifts::rangeLnProb( size_t i )
{
    double b = ranges[i].birth;
    double d = ranges[i].death;

    double present = times.front();

    size_t bi = findIndex(b);
    size_t di = findIndex(d);

    // include speciation density
    double lnProb = log( birth[bi] );

    // skip the rest for extant taxa with no fossil samples
    if ( taxa[i].getMaxAge() == present )
    {
        return lnProb;
    }

    // include extinction density
    if ( ranges[i].survived == false ) lnProb += log( death[di] );

    double psi_b_d = 0.0;

    for ( size_t j = di; j <= bi; j++ )
    {
        double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

        double dt = std::min(b, t_0) - std::max(d, times[j]);

        lnProb -= (birth[j] + death[j])*dt;

        psi_b_d += fossil[j]*dt;
    }

    // fossil non-sampling normalization over [d,b]
    lnProb -= psi_b_d;

    if ( condition == "sampling" )
    {
        lnProb -= log(-expm1(-psi_b_d));
    }

    return lnProb;
}
