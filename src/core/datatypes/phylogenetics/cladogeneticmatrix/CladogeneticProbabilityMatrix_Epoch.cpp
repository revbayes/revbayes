//
//  CladogeneticProbabilityMatrix_Epoch.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 10/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "CladogeneticProbabilityMatrix_Epoch.h"

using namespace RevBayesCore;

#include "RbException.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
CladogeneticProbabilityMatrix_Epoch::CladogeneticProbabilityMatrix_Epoch(void) :
CladogeneticProbabilityMatrix(0)
{
    ; // do nothing
}

CladogeneticProbabilityMatrix_Epoch::CladogeneticProbabilityMatrix_Epoch(size_t n) :
CladogeneticProbabilityMatrix( n )
{
    ; // do nothing
}

CladogeneticProbabilityMatrix_Epoch::~CladogeneticProbabilityMatrix_Epoch(void)
{
    ; // do nothing
}

CladogeneticProbabilityMatrix_Epoch* CladogeneticProbabilityMatrix_Epoch::clone( void ) const
{
    
    return new CladogeneticProbabilityMatrix_Epoch( *this );
}

size_t CladogeneticProbabilityMatrix_Epoch::findEpochIndex( double t ) const
{
    int i = 0;
    //    const std::vector<double> epoch_times = epochTimes->getValue();
    
    for (i = 0; i < epochTimes.size() && t < epochTimes[i]; i++)
    {
        ; // do nothing
    }
    return i;
}

const RbVector<double>& CladogeneticProbabilityMatrix_Epoch::getEpochTimes(void) const
{
    return epochTimes;
}

const RbVector<CladogeneticProbabilityMatrix>& CladogeneticProbabilityMatrix_Epoch::getCladogeneticProbabilityMatrix(void) const
{
    return epochCladogeneticProbabilityMatrices;
}

const CladogeneticProbabilityMatrix& CladogeneticProbabilityMatrix_Epoch::getCladogeneticProbabilityMatrix(double t) const
{
    size_t k = findEpochIndex(t);
    return epochCladogeneticProbabilityMatrices[k];
}

void CladogeneticProbabilityMatrix_Epoch::printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const {
    
    std::stringstream prev_time( "Inf" );
    for (size_t i = 0; i < epochTimes.size(); i++) {
        
        std::stringstream curr_time;
        curr_time << epochTimes[i];
        o << "\n";
        o << "Epoch with index " << i << "\n";
        o << "Epoch time from " << prev_time.str() << " to " << curr_time.str() << "\n";
        o << "Epoch cladogenetic probability matrix:\n";
        // append each clado matrix's print statement
        epochCladogeneticProbabilityMatrices[i].printForUser(o, sep, l, left);
        prev_time.str( curr_time.str() );
    }
}

json CladogeneticProbabilityMatrix_Epoch::toJSON() const
{
    json epochs;

    for (size_t i = 0; i < epochTimes.size(); i++)
    {
	json epoch;
	epoch["endTime"] = epochTimes[i];
	epoch["matrix"] = epochCladogeneticProbabilityMatrices[i].toJSON();
	epochs.push_back(epoch);
    }

    return epochs;
}

void CladogeneticProbabilityMatrix_Epoch::setEpochCladogeneticProbabilityMatrix(const RbVector<CladogeneticProbabilityMatrix>& cp)
{
    epochCladogeneticProbabilityMatrices = cp;
}

void CladogeneticProbabilityMatrix_Epoch::setEpochCladogeneticProbabilityMatrix(const CladogeneticProbabilityMatrix& cp, double t)
{
    size_t k = findEpochIndex(t);
    epochCladogeneticProbabilityMatrices[k] = cp;
}

void CladogeneticProbabilityMatrix_Epoch::setEpochTimes(const RbVector<double> &t)
{
    epochTimes = t;
}

std::map<std::vector<unsigned>, double> CladogeneticProbabilityMatrix_Epoch::getEventMap(double t)
{
    size_t k = findEpochIndex(t);
    return epochCladogeneticProbabilityMatrices[k].getEventMap();
}

const std::map<std::vector<unsigned>, double>& CladogeneticProbabilityMatrix_Epoch::getEventMap(double t) const
{
    size_t k = findEpochIndex(t);
    return epochCladogeneticProbabilityMatrices[k].getEventMap();
}

void CladogeneticProbabilityMatrix_Epoch::setEventMap(const std::map<std::vector<unsigned int>, double>& m, double t)
{
    size_t k = findEpochIndex(t);
    epochCladogeneticProbabilityMatrices[k].setEventMap(m);
}

void CladogeneticProbabilityMatrix_Epoch::update(void)
{
    ; // do nothing?
}

void CladogeneticProbabilityMatrix_Epoch::initFromString( const std::string &s )
{
    throw RbException("MJL (23/10/2016): Missing implementation!");
}

