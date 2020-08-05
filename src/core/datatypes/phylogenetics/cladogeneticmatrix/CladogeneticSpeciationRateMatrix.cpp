//
//  CladogeneticSpeciationRateMatrix.h
//
//  Created by Will Freyman on 8/1/17
//

#include "CladogeneticSpeciationRateMatrix.h"

#include <fstream>
#include <iomanip>
#include <utility>

#include "DagNode.h"
#include "RbException.h"

using namespace RevBayesCore;

CladogeneticSpeciationRateMatrix::CladogeneticSpeciationRateMatrix(void) :
num_states( 0 )
{
    ; // do nothing
}

CladogeneticSpeciationRateMatrix::CladogeneticSpeciationRateMatrix(size_t n) :
num_states( n )
{
    ; // do nothing
}

CladogeneticSpeciationRateMatrix::~CladogeneticSpeciationRateMatrix(void)
{
    ; // do nothing
}

CladogeneticSpeciationRateMatrix& CladogeneticSpeciationRateMatrix::assign(const Assignable &m)
{
    
    const CladogeneticSpeciationRateMatrix *cp = dynamic_cast<const CladogeneticSpeciationRateMatrix*>(&m);
    if ( cp != NULL )
    {
        return operator=(*cp);
    }
    else
    {
        throw RbException("Could not assign cladogenetic speciation rate matrix.");
    }
}


CladogeneticSpeciationRateMatrix* CladogeneticSpeciationRateMatrix::clone( void ) const
{
    
    return new CladogeneticSpeciationRateMatrix( *this );
}

void CladogeneticSpeciationRateMatrix::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, CladogeneticProbabilityMatrix &rv) const
{

    if ( n == "getCladogeneticProbabilityMatrix" )
    {
        rv = cladogenetic_probability_matrix;
    }
}

void CladogeneticSpeciationRateMatrix::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{

    if ( n == "getSpeciationRateSumPerState" )
    {
        rv = speciation_rate_sum_per_state;
    }
}

std::map<std::vector<unsigned>, double> CladogeneticSpeciationRateMatrix::getEventMap(double t)
{
    return event_map;
}

const std::map<std::vector<unsigned>, double>& CladogeneticSpeciationRateMatrix::getEventMap(double t) const
{
    return event_map;
}

std::vector<double> CladogeneticSpeciationRateMatrix::getSpeciationRateSumPerState(void)
{
    return speciation_rate_sum_per_state;
}

const std::vector<double>& CladogeneticSpeciationRateMatrix::getSpeciationRateSumPerState(void) const
{
    return speciation_rate_sum_per_state;
}

void CladogeneticSpeciationRateMatrix::setSpeciationRateSumPerState(std::vector<double> r)
{
    speciation_rate_sum_per_state = r;
}

CladogeneticProbabilityMatrix CladogeneticSpeciationRateMatrix::getCladogeneticProbabilityMatrix(void)
{
    return cladogenetic_probability_matrix;
}
const CladogeneticProbabilityMatrix& CladogeneticSpeciationRateMatrix::getCladogeneticProbabilityMatrix(void) const
{
    return cladogenetic_probability_matrix;
}

void CladogeneticSpeciationRateMatrix::setCladogeneticProbabilityMatrix(CladogeneticProbabilityMatrix p)
{
    cladogenetic_probability_matrix = p;
}

size_t CladogeneticSpeciationRateMatrix::getNumberOfStates( void ) const
{
    return num_states;
}

void CladogeneticSpeciationRateMatrix::initFromString( const std::string &s )
{
    throw RbException("WAF (8/1/2017): Missing implementation!");
}

void CladogeneticSpeciationRateMatrix::setEventMap(std::map<std::vector<unsigned>, double> m)
{
    event_map = m;
}


size_t CladogeneticSpeciationRateMatrix::size( void ) const
{
    return num_states;
}


void CladogeneticSpeciationRateMatrix::printForUser(std::ostream &o, const std::string &sep, int l, bool left) const
{
    std::streamsize previous_precision = o.precision();
    std::ios_base::fmtflags previous_flags = o.flags();
    
    o << "[ ";
    o << std::fixed;
    o << std::setprecision(4);

    o << "\n";
    std::map<std::vector<unsigned>, double>::const_iterator it;
    for (it = event_map.begin(); it != event_map.end(); it++)
    {
        if (it != event_map.begin())
            o << ",\n";
        const std::vector<unsigned>& v = it->first;
        o << "  " << "( " << v[0] << " -> " << v[1] << ", " << v[2] << " ) = " << it->second;
    }
    o << "\n]\n";
    
    o.setf(previous_flags);
    o.precision(previous_precision);
}




void CladogeneticSpeciationRateMatrix::printForSimpleStoring(std::ostream &o, const std::string &sep, int l, bool left, bool flatten) const
{
    /*
    
    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=0; i < size(); i++)
    {
        if (i > 0)
        {
            o << sep;
        }
        for (size_t j = 0; j < size(); ++j)
        {
            if (j > 0)
            {
                o << sep;
            }
            o << getRate( i, j, 1e-6, 1.0);
        }
        
    }
    */
}



void CladogeneticSpeciationRateMatrix::printForComplexStoring(std::ostream &o, const std::string &sep, int l, bool left, bool flatten) const
{
    /*
    o << "[ ";
    
    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=0; i < size(); i++)
    {
        o << "[ ";
        for (size_t j = 0; j < size(); ++j)
        {
            if (j != 0)
            {
                o << ", ";
            }
            o << getRate( i, j, 1e-6, 1.0);
        }
        o <<  " ]";
        
        if (i == size()-1)
        {
            o << " ]";
        }
        else
        {
            o << " ,";
        }
        
    }
     */
    
}
