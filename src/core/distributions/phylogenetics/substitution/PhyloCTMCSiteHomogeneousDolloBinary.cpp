#include "PhyloCTMCSiteHomogeneousDolloBinary.h"
#include "RbMathFunctions.h"

RevBayesCore::PhyloCTMCSiteHomogeneousDolloBinary::PhyloCTMCSiteHomogeneousDolloBinary(const TypedDagNode<Tree> *t, bool c, size_t nSites, bool amb, DolloAscertainmentBias::Coding ty, bool norm) :
    PhyloCTMCSiteHomogeneousDollo(  t, 1, c, nSites, amb, DolloAscertainmentBias::Coding(ty))
{

}

bool RevBayesCore::PhyloCTMCSiteHomogeneousDolloBinary::isSitePatternCompatible( std::map<size_t, size_t> charCounts )
{
    size_t zero = charCounts[0];
    size_t one  = charCounts[1];

    bool compatible = true;

    if (zero == 0 && (coding & DolloBinaryAscertainmentBias::NOPRESENCESITES) )
    {
        compatible = false;
    }
    else if (one == 0  && (coding & DolloBinaryAscertainmentBias::NOABSENCESITES) )
    {
        compatible = false;
    }
    else if (zero == 1 && one > 0 && (coding & DolloBinaryAscertainmentBias::NOSINGLETONABSENCE) )
    {
        compatible = false;
    }
    else if (zero > 0 && one == 1 && (coding & DolloBinaryAscertainmentBias::NOSINGLETONPRESENCE) )
    {
        compatible = false;
    }

    return compatible;
}

double RevBayesCore::PhyloCTMCSiteHomogeneousDolloBinary::computeIntegratedNodeCorrection(const std::vector<std::vector<std::vector<double> > >& partials, size_t node_index, size_t mask, size_t mixture, const std::vector<double> &f)
{
    /*
     * partials[a][b][c] contains the probability of observing:
     * autapomorphy pattern c
     * on background state b
     * conditioned on node state a
     */

    double prob = 0.0;

    if (coding & DolloBinaryAscertainmentBias::NOABSENCESITES)
        prob += partials[0][1][0];

    if (coding & DolloBinaryAscertainmentBias::NOSINGLETONPRESENCE)
        prob += partials[0][1][1];

    // if there is only one observed tip, then don't double-count singleton gains
    if ((coding & DolloBinaryAscertainmentBias::NOPRESENCESITES)
            && maskNodeObservationCounts[mask][node_index] == maskObservationCounts[mask]
            && (maskObservationCounts[mask] > 1 || !(coding & DolloBinaryAscertainmentBias::NOSINGLETONPRESENCE)) )
        prob += partials[0][0][0];

    // if there are only two observed tips, then don't double-count singleton gains
    // if there is only one observed tip, then don't double-count absence sites
    if ((coding & DolloBinaryAscertainmentBias::NOSINGLETONABSENCE)
            && (maskObservationCounts[mask] > 1 || !(coding & DolloBinaryAscertainmentBias::NOABSENCESITES))
            && (maskObservationCounts[mask] > 2 || !(coding & DolloBinaryAscertainmentBias::NOSINGLETONPRESENCE)) )
    {
        if (maskNodeObservationCounts[mask][node_index] == maskObservationCounts[mask])
            prob += partials[0][0][1];
        else if (maskNodeObservationCounts[mask][node_index] == maskObservationCounts[mask] - 1)
            prob += partials[0][0][0];
    }

    // impose a per-mixture boundary
    if (prob < 0.0 || prob > 1.0)
    {
        prob = RbConstants::Double::nan;
    }

    if (prob == 0.0 && !this->tau->getValue().getNode(node_index).isTip())
    {
        //prob = RbConstants::Double::nan;
    }

    prob = integrationFactors[mixture]*(1.0 - prob);

    return prob;
}
