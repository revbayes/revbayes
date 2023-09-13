//
//  RateMatrix_Biogeography.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/16/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

//#define DEBUG_DEC

#include <cstddef>
#include <cmath>
#include <list>
#include <complex>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_Biogeography.h"
#include "RbException.h"
#include "RbVector.h"
#include "TransitionProbabilityMatrix.h"
#include "Assignable.h"
#include "GeneralRateMatrix.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_Biogeography::RateMatrix_Biogeography(size_t ns, size_t nc, size_t mrs) : GeneralRateMatrix( ns ),
    numCharacters(nc),
    num_states(ns),
    useSquaring(ns > 32),
    hasEigenSystem(false),
    theEigenSystem(nullptr),
    dispersalRates( RbVector<RbVector<double > >( numCharacters, RbVector<double>(numCharacters, 1.0) ) ),
    extirpationRates( RbVector<double>( numCharacters, 1.0) ),
    scalingFactor(1.0),
    birthRate(0.0),
    rescaleMatrix(false),
    maxRangeSize(mrs),
    stationaryMatrix( TransitionProbabilityMatrix(num_states) ),
    accessedTransitionProbabilities( std::list<double>() ),
    maxSizeStoredTransitionProbabilites(1e3),
    useStoredTransitionProbabilities(true)
{
    
    makeBits();
    makeTransitions();
    
    for (size_t i = 0; i < num_states; ++i)
    {
        for (size_t j = 0; j < num_states; ++j)
        {
            (*the_rate_matrix)[i][j] = 0.0;
        }
    }
    update();

}


/** Copy constructor */
RateMatrix_Biogeography::RateMatrix_Biogeography(const RateMatrix_Biogeography& m) : GeneralRateMatrix( m ),
    stationaryMatrix(m.stationaryMatrix)
{
    
    bits                 = m.bits;
    inverseBits          = m.inverseBits;
    statesToBitsByNumOn  = m.statesToBitsByNumOn;
    bitsToStatesByNumOn  = m.bitsToStatesByNumOn;
    transitions          = m.transitions;
    lossOrGain           = m.lossOrGain;
    transitionAreas      = m.transitionAreas;
    numCharacters        = m.numCharacters;
    num_states           = m.num_states;
    useSquaring          = m.useSquaring;
    hasEigenSystem       = false;
    theEigenSystem       = nullptr;
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    dispersalRates       = m.dispersalRates;
    extirpationRates     = m.extirpationRates;
    useCladogenesis      = m.useCladogenesis;
    maxRangeSize         = m.maxRangeSize;
    rescaleMatrix        = m.rescaleMatrix;
    scalingFactor        = m.scalingFactor;
    storedTransitionProbabilities = m.storedTransitionProbabilities;
    accessedTransitionProbabilities = m.accessedTransitionProbabilities;
    maxSizeStoredTransitionProbabilites = m.maxSizeStoredTransitionProbabilites;
    useStoredTransitionProbabilities = m.useStoredTransitionProbabilities;
    changedAreas = m.changedAreas;
    affectingAreas = m.affectingAreas;
    
    update();

}


/** Destructor */
RateMatrix_Biogeography::~RateMatrix_Biogeography(void) {
    delete theEigenSystem;
}


RateMatrix_Biogeography& RateMatrix_Biogeography::operator=(const RateMatrix_Biogeography &r) {
    
    if (this != &r)
    {
        GeneralRateMatrix::operator=( r );
        
        delete theEigenSystem;
        
        bits                 = r.bits;
        inverseBits          = r.inverseBits;
        statesToBitsByNumOn  = r.statesToBitsByNumOn;
        bitsToStatesByNumOn  = r.bitsToStatesByNumOn;
        transitions          = r.transitions;
        lossOrGain           = r.lossOrGain;
        transitionAreas      = r.transitionAreas;
        numCharacters        = r.numCharacters;
        num_states           = r.num_states;
        useSquaring          = r.useSquaring;
        hasEigenSystem       = false;
        theEigenSystem       = nullptr;
        c_ijk                = r.c_ijk;
        cc_ijk               = r.cc_ijk;
        dispersalRates       = r.dispersalRates;
        extirpationRates     = r.extirpationRates;
        useCladogenesis      = r.useCladogenesis;
        maxRangeSize         = r.maxRangeSize;
        rescaleMatrix        = r.rescaleMatrix;
        scalingFactor        = r.scalingFactor;
        storedTransitionProbabilities = r.storedTransitionProbabilities;
        accessedTransitionProbabilities = r.accessedTransitionProbabilities;
        maxSizeStoredTransitionProbabilites = r.maxSizeStoredTransitionProbabilites;
        useStoredTransitionProbabilities = r.useStoredTransitionProbabilities;
        changedAreas = r.changedAreas;
        affectingAreas = r.affectingAreas;

        update();
        
    }
    
    return *this;
}

RateMatrix_Biogeography& RateMatrix_Biogeography::assign(const Assignable &m)
{
    const RateMatrix_Biogeography *rm = dynamic_cast<const RateMatrix_Biogeography*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }

}

double RateMatrix_Biogeography::averageRate(void) const
{
    double ave = 0.0;
    for (size_t i=0; i<num_states; i++)
        ave -= (*the_rate_matrix)[i][i];
    return ave / num_states;
}



void RateMatrix_Biogeography::fillRateMatrix( void )
{
        
    MatrixReal& m = *the_rate_matrix;
    
    std::vector<std::vector<unsigned> >::iterator it;
    std::vector<unsigned>::iterator jt;
    
    // i is the integer-valued index of the starting state
    for (size_t i = 0; i < transitions.size(); i++)
    {
        unsigned startState = (unsigned)i;
        
        // get range size weights
        int n = 0;
        for (size_t j = 0; j < statesToBitsByNumOn[i].size(); j++)
            n += statesToBitsByNumOn[i][j];
        
        double sum = 0.0;
        
        // j indexes of the possible moves leaving state i
        // transitions[i][j] gives the value of destination states j for i->j
        for (size_t j = 0; j < transitions[i].size(); j++)
        {
            unsigned endState = transitions[i][j];
            
            double v = 0.0;
            
            std::vector<unsigned> b1 = statesToBitsByNumOn[startState];
            std::vector<unsigned> b2 = statesToBitsByNumOn[endState];
            
//            std::cout << getRangeStr(b1) << "->" << getRangeStr(b2) << " " << (lossOrGain[i][j]==0 ? "LOSS" : "GAIN") << " :\n";
            
            // extinction
            if (lossOrGain[i][j] == 0)
            {
                unsigned changed_area = changedAreas[i][j];
                std::vector<unsigned> affecting_areas = affectingAreas[i][j];
                
                for (size_t k = 0; k < affecting_areas.size(); k++)
                {
                    double vt = extirpationRates[changed_area];
//                        std::cout << "\t" << vt << "dr[ " <<  affecting_areas[k] << " ][ " << changed_area << "]\n";
                    v += vt;
                }
            }
            // dispersal
            else if (lossOrGain[i][j] == 1)
            {
                unsigned changed_area = changedAreas[i][j];
                std::vector<unsigned> affecting_areas = affectingAreas[i][j];
                
                for (size_t k = 0; k < affecting_areas.size(); k++)
                {
                    double vt = dispersalRates[ affecting_areas[k] ][changed_area];
//                    std::cout << "\t" << vt << "dr[ " <<  affecting_areas[k] << " ][ " << changed_area << "]\n";
                    v += vt;
                }
            }
//            std::cout << "\tsum = " << sum << "\n";
                        
            // store value
            m[ startState ][ endState ] = v;
            
            // accumulate diagonal sum
            sum += v;

        }
        
        // set diagonal
        m[startState][startState] = -sum;
    }
//    std::cout << m << "\n";
    
    // set flags
    needs_update = true;
}

/** Do precalculations on eigenvectors */
void RateMatrix_Biogeography::calculateCijk(void)
{
    
    if ( theEigenSystem->isComplex() == false )
    {
        // real case
        const MatrixReal& ev  = theEigenSystem->getEigenvectors();
        const MatrixReal& iev = theEigenSystem->getInverseEigenvectors();
        double* pc = &c_ijk[0];
        for (size_t i=0; i<num_states; i++)
        {
            for (size_t j=0; j<num_states; j++)
            {
                for (size_t k=0; k<num_states; k++)
                {
                    *(pc++) = ev[i][k] * iev[k][j];
                }
            }
        }
    }
    else
    {
        // complex case
        const MatrixComplex& cev  = theEigenSystem->getComplexEigenvectors();
        const MatrixComplex& ciev = theEigenSystem->getComplexInverseEigenvectors();
        std::complex<double>* pc = &cc_ijk[0];
        for (size_t i=0; i<num_states; i++)
        {
            for (size_t j=0; j<num_states; j++)
            {
                for (size_t k=0; k<num_states; k++)
                {
                    *(pc++) = cev[i][k] * ciev[k][j];
                }
            }
        }
    }
}


/** Calculate the transition probabilities */
void RateMatrix_Biogeography::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    double t = scalingFactor * rate * (startAge - endAge);
    
    if (t != 0.0) {
        double digits = 8;
        double factor = std::pow(10.0, digits - std::ceil(std::log10(std::fabs(t))));

        t = round(t * factor) / factor;
    }
    
    
    // Do we already have P(t)?
    std::map<double, TransitionProbabilityMatrix>::const_iterator it = storedTransitionProbabilities.find(t);
    bool found = it != storedTransitionProbabilities.end();
    if (found) {
        
        // update the transition probs
        P = it->second;
        
        // this time was most recently accessed
        accessedTransitionProbabilities.remove(it->first);
        accessedTransitionProbabilities.push_front(it->first);
    }
    else {
        
        if ( useSquaring || true )
        {
            //We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
            exponentiateMatrixByScalingAndSquaring(t, P);
        }
        else if ( theEigenSystem->isComplex() == false )
        {
            tiProbsEigens(t, P);
        }
        else
        {
            tiProbsComplexEigens(t, P);
        }

        if (useStoredTransitionProbabilities) {
            storedTransitionProbabilities.insert( std::pair<double, TransitionProbabilityMatrix>(t, P) );
            accessedTransitionProbabilities.push_front(t);
        }
        if (accessedTransitionProbabilities.size() > maxSizeStoredTransitionProbabilites)
        {
            double rem_time = accessedTransitionProbabilities.back();
            storedTransitionProbabilities.erase(rem_time);
            accessedTransitionProbabilities.pop_back();
        }
//        std::cout << storedTransitionProbabilities.size() << "\n";
    }
    
    return;
}

RateMatrix_Biogeography* RateMatrix_Biogeography::clone( void ) const
{
    return new RateMatrix_Biogeography( *this );
}

const RbVector<RbVector<double> >& RateMatrix_Biogeography::getDispersalRates(void) const
{
    return dispersalRates;
}

const RbVector<double>& RateMatrix_Biogeography::getExtirpationRates(void) const
{
    return extirpationRates;
}

std::string RateMatrix_Biogeography::getRangeStr(const std::vector<unsigned>& v)
{
    std::stringstream ss;
    for (size_t j = 0; j < v.size(); j++)
        ss << v[j];
    return ss.str();
}

std::vector<double> RateMatrix_Biogeography::getStationaryFrequencies(void) const
{
    // MJL: The initial DEC model uses flat stationary frequencies.
    //      Two alternative solutions in RevBayes (not yet exposed in RevLanguage):
    //          1) DEC conditioned on survival without cladogenesis
    //          2) DEC conditioned on survival with cladogensis
    std::vector<double> f(num_states, 1.0/num_states);
    return(f);
}

void RateMatrix_Biogeography::makeBits(void)
{
    
    size_t num_all_states = (size_t)pow(2,numCharacters);
    bitsByNumOn.resize(numCharacters+1);
    
    
//    if (includeNullRangeAsState == true) {
//        statesToBitsByNumOn.resize(num_all_states);
//    }
//    else {
        statesToBitsByNumOn.resize(num_all_states-1);
//    }
    bits = std::vector<std::vector<unsigned> >(num_all_states, std::vector<unsigned>(numCharacters, 0));
    
//    if (includeNullRangeAsState) {
//        bitsByNumOn[0].push_back(bits[0]);
//    }
    for (size_t i = 1; i < num_all_states; i++)
    {
        size_t m = i;
        for (size_t j = 0; j < numCharacters; j++)
        {
            bits[i][j] = m % 2;
            m /= 2;
            if (m == 0)
                break;
        }
        size_t j = numBitsOn(bits[i]);
        bitsByNumOn[j].push_back(bits[i]);
        
    }
    for (size_t i = 0; i < num_all_states; i++)
    {
        inverseBits[ bits[i] ] = (unsigned)i;
    }
    
    // assign state to each bit vector, sorted by numOn
    size_t k = 0;
    for (size_t i = 0; i < bitsByNumOn.size(); i++)
    {
        for (size_t j = 0; j < bitsByNumOn[i].size(); j++)
        {
            statesToBitsByNumOn[k++] = bitsByNumOn[i][j];
        }
    }
    
    for (size_t i = 0; i < statesToBitsByNumOn.size(); i++)
    {
        bitsToStatesByNumOn[ statesToBitsByNumOn[i] ] = (unsigned)i;
    }
    

    
}

size_t RateMatrix_Biogeography::numBitsOn(std::vector<unsigned> v)
{
    size_t n = 0;
    for (size_t i = 0; i < v.size(); i++) {
        n += v[i];
    }
    return n;
}

void RateMatrix_Biogeography::makeTransitions(void)
{
    
    transitions.resize(num_states);
    lossOrGain.resize(num_states);
    transitionAreas.resize(num_states);
    changedAreas.resize(num_states);
    affectingAreas.resize(num_states);
    
    // populate integer-valued transitions between states
    for (size_t i = 0; i < num_states; i++)
    {
        std::vector<unsigned> b = statesToBitsByNumOn[i];
        
        // each row has b.size() events (excluding i==0)
        for (size_t j = 0; j < b.size(); j++)
        {
            std::vector<unsigned> tmp = b;
            
            // change the range cfg at area j
            tmp[j] = (b[j] == 0 ? 1 : 0);
            
            // ignore events larger than maxRangeSize
            if (numBitsOn(tmp) > maxRangeSize || numBitsOn(b) > maxRangeSize)
            {
                continue;
            }
            // ignore all-0 ranges (this is an extinction event)
            if (numBitsOn(tmp) == 0) {
                continue;
            }
            
            // store integer-valued event
            transitions[i].push_back(bitsToStatesByNumOn[tmp]);
            
            // is event a gain or a loss?
            lossOrGain[i].push_back(tmp[j]);
            
            std::vector<unsigned> a;

            changedAreas[i].push_back((unsigned)j);
            for (size_t k = 0; k < b.size(); k++)
            {
                if (b[k]==1)
                {
                    a.push_back((unsigned)k);
                }
            }
            affectingAreas[i].push_back(a);
        }
    }
}

void RateMatrix_Biogeography::setDispersalRates(const RbVector<RbVector<double> >& dr)
{
    dispersalRates = dr;
    needs_update = true;
}

void RateMatrix_Biogeography::setExtirpationRates(const RbVector<double>& er)
{
    extirpationRates = er;
    needs_update = true;
}

/** Calculate the transition probabilities for the real case */
void RateMatrix_Biogeography::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = theEigenSystem->getRealEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<double> eigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
    {
        eigValExp[s] = exp(eigenValue[s] * t);
    }
    
    // calculate the transition probabilities
    const double* ptr = &c_ijk[0];
    double*         p = P.theMatrix;
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++, ++p)
        {
            double sum = 0.0;
            for (size_t s=0; s<num_states; s++)
            {
                sum += (*ptr++) * eigValExp[s];
            }
            (*p) = (sum < 0.0) ? 0.0 : sum;
        }
    }
}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_Biogeography::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const {
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValueReal = theEigenSystem->getRealEigenvalues();
    const std::vector<double>& eigenValueComp = theEigenSystem->getImagEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<std::complex<double> > ceigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
    {
        std::complex<double> ev = std::complex<double>(eigenValueReal[s], eigenValueComp[s]);
        ceigValExp[s] = exp(ev * t);
    }
    
    // calculate the transition probabilities
    const std::complex<double>* ptr = &cc_ijk[0];
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
            std::complex<double> sum = std::complex<double>(0.0, 0.0);
            for (size_t s=0; s<num_states; s++)
                sum += (*ptr++) * ceigValExp[s];
            P[i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
        }
    }
}

/** Create the eigen system */
void RateMatrix_Biogeography::initializeEigenSystem(void) {
    
    if (!hasEigenSystem) {
        
        // make the eigen system
        theEigenSystem = new EigenSystem(the_rate_matrix);
        c_ijk.resize(num_states * num_states * num_states);
        cc_ijk.resize(num_states * num_states * num_states);

        // flag as created
        hasEigenSystem = true;
        
    }

}


/** Update the eigen system */
void RateMatrix_Biogeography::updateEigenSystem(void) {
    
    // initialize, if hasn't been done already
    initializeEigenSystem();

    // update the eigen system
    theEigenSystem->update();
    calculateCijk();
    
}


void RateMatrix_Biogeography::update( void ) {
    
    if ( needs_update )
    {
        // assign all rate matrix elements
        fillRateMatrix();
        
        // rescale
        scalingFactor = 1.0;
        if (rescaleMatrix)
            rescaleToAverageRate(1.0);
        
        if (!useSquaring)
            // get transition probs
            updateEigenSystem();
        
        // clear the stored transition probabilities
//        std::cout << "Clearing " << storedTransitionProbabilities.size() << "\n";
        storedTransitionProbabilities.clear();
        accessedTransitionProbabilities = std::list<double>();
        
        // clean flags
        needs_update = false;
    }
}

