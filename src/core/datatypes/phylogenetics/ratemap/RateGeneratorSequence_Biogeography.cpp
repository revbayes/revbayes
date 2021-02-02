#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "BranchHistory.h"
#include "CharacterEventDiscrete.h"
#include "RateGeneratorSequence_Biogeography.h"
#include "Cloneable.h"
#include "DistanceRateModifier.h"
#include "RateGenerator.h"
#include "RateGeneratorSequence.h"
#include "RateGeneratorSequenceUsingMatrix.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;

RateGeneratorSequence_Biogeography::RateGeneratorSequence_Biogeography(size_t nc, bool fe, unsigned mrs) : RateGeneratorSequenceUsingMatrix(2, nc)
{

    forbidExtinction = fe;
    distanceRateModifier = NULL;
    distancePower = 0.0;
    maxRangeSize = mrs;
    if (mrs == 0 || mrs > nc)
        mrs = (unsigned)nc;
//    rootFrequencies = std::vector<double>(2,0.5);


    epochs = std::vector<double>(1,0.0);
    numEpochs = epochs.size();
    epochOffset = numEpochs * num_characters;
//    extinctionValues = std::vector<double>(numEpochs * this->numCharacters, 1.0);
//    dispersalValues = std::vector<double>(numEpochs * this->numCharacters * this->numCharacters, 1.0);
    adjacentAreaVector = std::vector<double>(numEpochs * this->num_characters * this->num_characters, 1.0);
    availableAreaVector = std::vector<double>(numEpochs * this->num_characters, 1.0);

    useDistanceRateModifier = false;
    useAreaAvailable = false;
    useAreaAdjacency = false;
    useDistanceDependence = false;
    useUnnormalizedRates = false;
//    useRootFrequencies = true;

    branchOffset=1;
    epochOffset=1;

}

RateGeneratorSequence_Biogeography::RateGeneratorSequence_Biogeography(const RateGeneratorSequence_Biogeography& m) : RateGeneratorSequenceUsingMatrix( m )
{
    distancePower = m.distancePower;
    maxRangeSize = m.maxRangeSize;

    epochs = m.epochs;
    numEpochs = m.numEpochs;
    epochOffset = m.epochOffset;
    adjacentAreaVector = m.adjacentAreaVector;
    availableAreaVector = m.availableAreaVector;

    useAreaAvailable = m.useAreaAvailable;
    useAreaAdjacency = m.useAreaAdjacency;
    useDistanceDependence = m.useDistanceDependence;
    useRootFrequencies = m.useRootFrequencies;
    useUnnormalizedRates = m.useUnnormalizedRates;

    distanceRateModifier     = m.distanceRateModifier;
    useDistanceRateModifier = m.useDistanceRateModifier;


    forbidExtinction = m.forbidExtinction;

    branchOffset = m.branchOffset;
    epochOffset = m.epochOffset;
}

RateGeneratorSequence_Biogeography::~RateGeneratorSequence_Biogeography(void)
{

}

RateGeneratorSequence_Biogeography& RateGeneratorSequence_Biogeography::operator=(const RateGeneratorSequence_Biogeography& r)
{
    if (this != &r)
    {
        RateGeneratorSequence::operator=( r );

        distancePower = r.distancePower;
        maxRangeSize = r.maxRangeSize;

        epochs = r.epochs;
        numEpochs = r.numEpochs;
        epochOffset = r.epochOffset;
        adjacentAreaVector = r.adjacentAreaVector;
        availableAreaVector = r.availableAreaVector;

        distanceRateModifier = r.distanceRateModifier;
        useDistanceRateModifier = r.useDistanceRateModifier;

        useAreaAdjacency = r.useAreaAdjacency;
        useAreaAvailable = r.useAreaAvailable;
        useDistanceDependence = r.useDistanceDependence;
        useUnnormalizedRates = r.useUnnormalizedRates;

        forbidExtinction = r.forbidExtinction;

        branchOffset = r.branchOffset;
        epochOffset = r.epochOffset;
    }
    return *this;
}

void RateGeneratorSequence_Biogeography::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{

    double branchLength = startAge - endAge;

    const RateGenerator* rm = rateMatrix;

//  @Michael: Please look at this again!!!
//    if (node.isRoot())
//        branchLength = node.getAge() * 5;

    double expPart = exp( -( rm->getRate(1,0, startAge, rate) + rm->getRate(0,1, startAge, rate) ) * branchLength);
    double p = rm->getRate(1,0, startAge, rate) / (rm->getRate(1,0, startAge, rate) + rm->getRate(0,1, startAge, rate));
    double q = 1.0 - p;

    P[0][0] = p + q * expPart;
    P[0][1] = q - q * expPart;
    P[1][0] = p - p * expPart;
    P[1][1] = q + p * expPart;
}


// @MJL: Fix this later 02/01/17
/*

void RateGeneratorSequence_Biogeography::calculateTransitionProbabilities(TransitionProbabilityMatrix& P, double startAge, double endAge, double r, size_t charIdx) const
{

    double currAge = startAge;

    const RateGenerator* rm = rateMatrix;

    // start at earliest epoch
    int epochIdx = getEpochIndex(startAge);


    // initialize P = I
    P[0][0] = 1.0;
    P[0][1] = 0.0;
    P[1][0] = 0.0;
    P[1][1] = 1.0;



    bool stop = false;
    while (!stop)
    {
        // get dispersal and extinction rates for site
        size_t idx = this->numCharacters * epochIdx + charIdx;

        double dispersalRate  = ( (availableAreaVector[ idx ] > 0) ? 1.0 : 0.0);
        double extinctionRate = ( (availableAreaVector[ idx ] > 0) ? 1.0 : 1e10);

        // get age of start of next oldest epoch
        double epochAge = epochs[ epochIdx ];

        // get next oldest age boundary (node or epoch)
        double incrAge = epochAge;

        // no more epochs, exit loop after computation
        if (endAge >= epochAge)
        {
            incrAge = endAge;
            stop = true;
        }

        // get branch length
        double diffAge = currAge - incrAge;

        // transition probabilities w/ sum-product
        double glr0 = rm->getRate(1,0, startAge, 1.0);
        double glr1 = rm->getRate(0,1, startAge, 1.0);
//        double glr0 = rm->getRate(1,0,currAge,extinctionRate);
//        double glr1 = rm->getRate(0,1,currAge,dispersalRate);
        double expPart = exp( -(glr0 + glr1) * r * diffAge);
        double p = glr0 / (glr0 + glr1);
        double q = 1.0 - p;

        double p00 = (p + q * expPart);
        double p01 = (q - q * expPart);
        double p10 = (p - p * expPart);
        double p11 = (q + p * expPart);

        double q00 = P[0][0];
        double q01 = P[0][1];
        double q10 = P[1][0];
        double q11 = P[1][1];

        P[0][0] = p00 * q00 + p01 * q10;
        P[0][1] = p00 * q01 + p01 * q11;
        P[1][0] = p10 * q00 + p11 * q10;
        P[1][1] = p10 * q01 + p11 * q11;

//        std::cout << P[0][0] << "," << P[0][1] << ";" << P[1][0] << "," << P[1][1] << "\n";



        if (!stop)
        {
            epochIdx += 1;
            currAge = epochAge;
        }
    }
}
 */

RateGeneratorSequence_Biogeography* RateGeneratorSequence_Biogeography::clone(void) const
{
    return new RateGeneratorSequence_Biogeography( *this );
}


//void RateGeneratorSequence_Biogeography::calculateTransitionProbabilities(TransitionProbabilityMatrix &P, double age) const
//{
//    const RateGenerator* rm = rateMatrix;
//
//    rm->calculateTransitionProbabilities(age, 0, 1.0, P);
//}


double RateGeneratorSequence_Biogeography::getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, unsigned* count, double age, double r) const
{
    size_t s = to->getState();

    if (static_cast<CharacterEventDiscrete*>(from[ to->getSiteIndex() ])->getState() == to->getState())
    {
//        std::cout << count[0] << " " << count[1] << "\n";
//        std::cout << node.getIndex() << " problem...\n";
        ;
    }

    // rate to extinction cfg is 0
    if (count[1] == 1 && s == 0 && forbidExtinction)
        return 0.0;

    // @Michael: Please look at this again too!!
    // rate for dispersal to exceed max range size is 0
//    if (!node.isRoot() && count[1] >= maxRangeSize && s == 1)
    if ( count[1] >= maxRangeSize && s == 1)
    {

        return 0.0;
    }

    // rate according to binary rate matrix Q(node)
//    double rate = rateMatrix->getRate(!s,s,age,r);
    double rate = rateMatrix->getRate(!s,s,age,r);

    // apply rate modifiers
    if (useDistanceRateModifier) // want this to take in age as an argument...
        rate *= distanceRateModifier->computeRateMultiplier(from, to, age);


    return rate;

}

double RateGeneratorSequence_Biogeography::getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age, double r) const
{
    unsigned n1 = (unsigned)numOn(from);
    unsigned n0 = (unsigned)(num_characters - n1);
    unsigned counts[2] = { n0, n1 };
    return getRate(from, to, counts, age, r);
}




double RateGeneratorSequence_Biogeography::getRate(size_t from, size_t to, double age, double r) const
{

    throw RbException("Missing implementation in RateGeneratorSequence_Biogeography (Sebastian)");

    return 0.0;
}

double RateGeneratorSequence_Biogeography::getSiteRate(CharacterEventDiscrete* from, CharacterEventDiscrete* to, double age, double r) const
{
    double rate = 0.0;
    size_t s = to->getState();
//    int charIdx = to->getIndex();
//    int epochIdx = getEpochIndex(age);

    // rate according to binary rate matrix Q(node)
    rate = rateMatrix->getRate(!s,s,age,r);
//    rate = rateMatrix->getRate(!s,s,age,r);

    // area effects
    if (useDistanceRateModifier)
        rate *= distanceRateModifier->computeSiteRateMultiplier(from,to,age);


    return rate;
}

double RateGeneratorSequence_Biogeography::getSiteRate( size_t from, size_t to, size_t charIdx, double age, double r) const
{
    double rate = 0.0;
    size_t s = to;

//    int pres = !s;
//    int epochIdx = getEpochIndex(age);

    // rate according to binary rate matrix Q(node)
    rate = rateMatrix->getRate(!s,s,age,r);


    // area effects
    if (useDistanceRateModifier)
    {
        double rm = distanceRateModifier->computeSiteRateMultiplier(from,to,charIdx,age);

        rate *= rm;
    }
//        rate *= distanceRateModifier->computeSiteRateModifier(node,from,to,charIdx,age);


    return rate;
}

double RateGeneratorSequence_Biogeography::getSumOfRates( std::vector<CharacterEvent*> from, std::vector<size_t> counts, double age, double r) const
{

    if (useUnnormalizedRates)
        return getUnnormalizedSumOfRates( from, counts, age);

    // get rate away away from currState
    unsigned n0 = (unsigned)counts[0];
//    if (useDistanceRateModifier)
//        n0 = distanceRateModifier->getNumAvailableAreas(node, from, age);
    unsigned n1 = (unsigned)counts[1];

    // forbid extinction events
//    if (counts[1] == 1 && forbidExtinction)
//        n1 = 0;
    if (counts[1] == 0 && forbidExtinction)
    {
        return 0.0;
    }

    // do not allow dispersal events if range size is at max
//    if (!node.isRoot() && counts[1] >= maxRangeSize)
    if ( counts[1] >= maxRangeSize)
    {

        n0 = 0;
    }


    // get (effective) num characters in each state
    double r0 = n1;
    double r1 = n0;
    if (useAreaAdjacency)
    {
        r1 = n1;
//        r0 = distanceRateModifier->getNumAvailableAreas(node,from,age);
        r0 = distanceRateModifier->getNumEmigratableAreas(from,age);
    }

//    if (useAreaAdjacency || useAreaAvailable)
//        r0 = distanceRateModifier->getNumAvailableAreas(node,from,age);
//
    // apply ctmc for branch
    r0 *= rateMatrix->getRate(1,0,age,r);
    r1 *= rateMatrix->getRate(0,1,age,r);



//    if (useRootFrequencies)
//    {
//        // root freqs
//        int countDiff[2][2] = { {1, -1}, {-1, 1} };
//        double rfr0 = pow(2*rootFrequencies[0], counts[0] + countDiff[0][0]) * pow(2*rootFrequencies[1], counts[1] + countDiff[0][1]);
//        double rfr1 = pow(2*rootFrequencies[0], counts[0] + countDiff[1][0]) * pow(2*rootFrequencies[1], counts[1] + countDiff[1][1]);
////        double rfr0 = rootFrequencies[0];
////        double rfr1 = rootFrequencies[1];
//        r0 *= rfr0;
//        r1 *= rfr1;
//    }

    // apply rate for branch.
    double sum = r0 + r1;

    return sum;
}

double RateGeneratorSequence_Biogeography::getSumOfRates( std::vector<CharacterEvent*> from, double age, double r) const
{
    size_t n1 = numOn(from);

    std::vector<size_t> c;
    c.push_back(num_characters - n1);
    c.push_back(n1);

    return RateGeneratorSequence_Biogeography::getSumOfRates( from, c, r, age);
}

double RateGeneratorSequence_Biogeography::getUnnormalizedSumOfRates( std::vector<CharacterEvent*> from, std::vector<size_t> counts, double age, double r) const
{

    size_t epochIdx = getEpochIndex(age);

    // apply ctmc for branch
    const RateGenerator* rm = rateMatrix;

    // get sum of rates
    double sum = 0.0;
    for (size_t i = 0; i < from.size(); i++)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(from[i])->getState();
        double v = availableAreaVector[ epochIdx * this->num_characters + i ];

        if (forbidExtinction && s == 1 && counts[1] == 0)
            sum += 0.0;
        else if (s == 1 && v > 0)
            sum += rm->getRate(1,0,age,1.0);
        else if (s == 1 && v == 0)
            sum += 1e10;
        else if (s == 0)
            sum += rm->getRate(0,1,age,1.0) * v;
    }

    // apply rate for branch
    sum *= r;

    return sum;
}



void RateGeneratorSequence_Biogeography::updateMap(void)
{
    if (needsUpdate)
    {
        ; // do nothing ...

        needsUpdate = false;
    }
}

double RateGeneratorSequence_Biogeography::getDistancePower(void) const
{
    return distancePower;
}

void RateGeneratorSequence_Biogeography::setDistancePower(double d)
{
    distancePower = d;
}

void RateGeneratorSequence_Biogeography::setDistanceRateModifier(const DistanceRateModifier& gdrm)
{
    useDistanceRateModifier = true;

    // ugly hack, prob better way to handle constness...
    distanceRateModifier = &const_cast<DistanceRateModifier&>(gdrm);

    epochs = distanceRateModifier->getEpochs();
    numEpochs = epochs.size();
    adjacentAreaVector = distanceRateModifier->getAdjacentAreaVector();
    availableAreaVector = distanceRateModifier->getAvailableAreaVector();
    useAreaAdjacency = distanceRateModifier->getUseAreaAdjacency();
    useAreaAvailable = distanceRateModifier->getUseAreaAvailable();
}

void RateGeneratorSequence_Biogeography::setGeographicDistancePowers(const DistanceRateModifier& gdrm)
{
    useDistanceRateModifier = true;
    distanceRateModifier->setGeographicDistancePowers(gdrm.getGeographicDistancePowers());
}

const DistanceRateModifier& RateGeneratorSequence_Biogeography::getDistanceRateModifier(void)
{
    return *distanceRateModifier;
}

const bool RateGeneratorSequence_Biogeography::isAreaAvailable(size_t charIdx, double age) const
{
    size_t epochIdx = getEpochIndex(age);
    return availableAreaVector[epochIdx*this->num_characters + charIdx] > 0.0;
}

const bool RateGeneratorSequence_Biogeography::areAreasAdjacent(size_t fromCharIdx, size_t toCharIdx, double age) const
{
    size_t epochIdx = getEpochIndex(age);
    return adjacentAreaVector[epochIdx*epochOffset + this->num_characters*fromCharIdx + toCharIdx] > 0.0;

}

const std::set<size_t> RateGeneratorSequence_Biogeography::getRangeAndFrontierSet(BranchHistory* bh, double age) const
{
    std::set<size_t> ret;
    const std::vector<CharacterEvent*>& from = bh->getParentCharacters();
//    const std::vector<CharacterEvent*>& to = bh->getParentCharacters();
    size_t epochIdx = getEpochIndex(age);
    const std::vector<std::set<size_t> >& adjacentAreaSet = distanceRateModifier->getAdjacentAreaSet();

//    std::set<size_t>
    for (size_t i = 0; i < from.size(); i++)
    {
        if (static_cast<CharacterEventDiscrete*>(from[i])->getState() == 1)
        {
            ret.insert(i);
            const std::set<size_t> adj = adjacentAreaSet[epochIdx*from.size() + i];
            std::set<size_t>::const_iterator it_adj;
            for (it_adj = adj.begin(); it_adj != adj.end(); it_adj++)
            {
                if (static_cast<CharacterEventDiscrete*>(from[*it_adj])->getState() == 0)
                {
                    ret.insert(*it_adj);
                }
            }
        }
    }

    return ret;
}

const std::vector<double>& RateGeneratorSequence_Biogeography::getEpochs(void) const
{
    return epochs;
}

size_t RateGeneratorSequence_Biogeography::numOn(const std::vector<CharacterEvent*>& s) const
{
    size_t n = 0;
    for (size_t i = 0; i < s.size(); i++)
        if (static_cast<CharacterEventDiscrete*>(s[i])->getState() == 1)
            n++;
    return n;
}

unsigned RateGeneratorSequence_Biogeography::getEpochIndex(double age) const
{
    // epochs are ordered from oldest to youngest, typically over (-neginf, 0.0)
    unsigned index = 0;
    while (age <= epochs[index] && index < epochs.size() - 1)
    {
        index++;
    };
    return index;
}
