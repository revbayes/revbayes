#include <fstream>
#include <string>
#include <iomanip>
#include <cstddef>
#include <vector>

#include "CharacterEventDiscrete.h"
#include "RateGenerator.h"
#include "RateMatrix.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDagNode.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "Simplex.h"

namespace RevBayesCore { class CharacterEvent; }
namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateGenerator::RateGenerator(size_t n) :
    num_states( n )
{
    ; // do nothing
}

RateGenerator::~RateGenerator(void)
{
    ; // do nothing
}


void RateGenerator::calculateTransitionProbabilities(double t, TransitionProbabilityMatrix& P) const
{


    calculateTransitionProbabilities(t, 0.0, 1.0, P);
}

size_t RateGenerator::getNumberOfStates( void ) const
{
    return num_states;
}

RbVector<std::string> RateGenerator::getStateDescriptions( void ) const
{
    RbVector<std::string> state_descriptions;
    for (size_t i = 0; i < num_states; i++) {
        std::stringstream ss;
        ss << i;
        state_descriptions.push_back(ss.str());
    }
    return state_descriptions;
}

double RateGenerator::getSumOfRates(std::vector<CharacterEvent*> from, const std::vector<size_t> &counts, double age, double rate) const
{

    // get the rate of leaving the sequence-state
    double sum = 0.0;
    for (size_t i = 0; i < num_states; ++i)
    {
        sum += -getRate(i, i, age, rate) * counts[i];
    }

    return sum;
}

double RateGenerator::getSumOfRates(std::vector<CharacterEvent*> from, double age, double rate) const
{

    // need dynamic allocation
    std::vector<size_t> counts = std::vector<size_t>(num_states,0);

    for (size_t i = 0; i < from.size(); ++i)
    {
        ++counts[ static_cast<CharacterEventDiscrete*>(from[i])->getState() ];
    }

    return getSumOfRates(from, counts, age, rate);
}

double RateGenerator::getSumOfRatesDifferential(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age, double rate) const
{

    double r = 0.0;
    
    size_t old_state = static_cast<CharacterEventDiscrete*>(from[ to->getSiteIndex() ])->getState();
    size_t new_state = to->getState();

    
    for (size_t possible_state = 0; possible_state < num_states; possible_state++)
    {
        // subtract the contribution of rates leaving the old state
        if (possible_state != old_state)
        {
            r -= getRate(old_state, possible_state, age, rate);
        }
        
        // add the contribution of rates leaving the new state
        if (possible_state != new_state)
        {
            r += getRate(new_state, possible_state, age, rate);
        }
    }
    
    return r;
}

void RateGenerator::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<double> >& rv) const
{

    // clear old values
    rv.clear();

    TransitionProbabilityMatrix P(num_states);

    double rate = static_cast<const TypedDagNode<double> *>( args[0] )->getValue();
    double start_age = static_cast<const TypedDagNode<double> *>( args[1] )->getValue();
    double end_age = static_cast<const TypedDagNode<double> *>( args[2] )->getValue();

    if (start_age < end_age)
    {
        double temp = start_age;
        start_age = end_age;
        end_age = temp;
    }

    calculateTransitionProbabilities( start_age, end_age, rate, P);

    for (size_t i = 0; i < num_states; i++)
    {
        RbVector<double> v;
        for (size_t j =0; j < num_states; j++)
        {
            v.push_back(P[i][j]);
        }
        rv.push_back(v);
    }

}

//void RateGenerator::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<double> >& rv) const
//{
//
//    // clear old values
//    rv.clear();
//    
//    TransitionProbabilityMatrix P(num_states);
//    
//    double rate = static_cast<const TypedDagNode<double> *>( args[0] )->getValue();
//    double start_age = static_cast<const TypedDagNode<double> *>( args[1] )->getValue();
//    double end_age = static_cast<const TypedDagNode<double> *>( args[2] )->getValue();
//    
//    calculateTransitionProbabilities( start_age, end_age, rate, P);
//    
//    for (size_t i = 0; i < num_states; i++)
//    {
//        RbVector<double> v;
//        for (size_t j =0; j < num_states; j++)
//        {
//            v.push_back(P[i][j]);
//        }
//        rv.push_back(v);
//    }
//    
//}


void RateGenerator::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const
{
    
    if ( n == "[]")
    {
        size_t n_states = this->getNumberOfStates();
        //    rv.resize(n_states);
        rv.clear();
        
        size_t from_idx = static_cast<const TypedDagNode<long> *>( args[0] )->getValue()-1;
        
        for (size_t to_idx = 0; to_idx < n_states; to_idx++)
        {
            rv.push_back(this->getRate(from_idx, to_idx, 0.0, 1.0));
        }
    }

}


void RateGenerator::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, Simplex &rv) const
{
    
    if (n == "getStationaryFrequencies")
    {
        // clear old values
        rv.clear();
        
        const RateMatrix *rm = dynamic_cast<const RateMatrix*>( this );
        if ( rm != NULL )
        {
            rv = rm->getStationaryFrequencies();
        }
        else
        {
            throw RbException("Stationary frequencies of RateGenerators that are not RateMatrices are not clearly defined.");
        }
    }
    
}

//void RateGenerator::executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const
//{
//    size_t n_states = this->getNumberOfStates();
////    rv.resize(n_states);
//    rv.clear();
//
//    size_t from_idx = static_cast<const TypedDagNode<long> *>( args[0] )->getValue()-1;
//
//    for (size_t to_idx = 0; to_idx < n_states; to_idx++)
//    {
//        rv.push_back(this->getRate(from_idx, to_idx, 0.0, 1.0));
//    }
//}


size_t RateGenerator::size( void ) const
{
    return num_states;

}

bool RateGenerator::simulateStochasticMapping(double startAge, double endAge, double rate, std::vector<size_t>& transition_states, std::vector<double>& transition_times) const
{
    throw RbException("simulateStochasticMapping not defined for abstract RateGenerator objects");
    return false;
}



json RateGenerator::toJSON() const
{
    json matrix;
    
    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=0; i < size(); i++)
    {
	json row;
        for (size_t j = 0; j < size(); ++j)
	    row.push_back( getRate( i, j, 1e-6, 1.0) );
	matrix.push_back(row);
    }

    return matrix;
}

void RateGenerator::printForUser(std::ostream &o, const std::string &sep, int l, bool left) const
{
    std::streamsize previous_precision = o.precision();
    std::ios_base::fmtflags previous_flags = o.flags();
    
    o << std::fixed;
    o << std::setprecision(4);

    // print the RbMatrix with each column of equal width and each column centered on the decimal
    for (size_t i=0; i < size(); i++)
    {
        if (i == 0)
        {
            o << "[ ";
        }
        else
        {
            o << "  ";
        }
        
        o << "[ ";
        for (size_t j = 0; j < size(); ++j)
        {
            if (j != 0)
            {
                o << ", ";
            }
            o << getRate( i, j, 1e-6,1.0);
        }
        o <<  " ]";
        
        if (i == size()-1)
        {
            o << " ]";
        }
        else
        {
            o << " ,\n";
        }
    }
    
    o.setf(previous_flags);
    o.precision(previous_precision);
    
}



void RateGenerator::printForSimpleStoring(std::ostream &o, const std::string &sep, int l, bool left, bool flatten) const
{
    
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
}



void RateGenerator::printForComplexStoring(std::ostream &o, const std::string &sep, int l, bool left, bool flatten) const
{
    
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
    
}
