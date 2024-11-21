#include <cmath>
#include <cstdint>
#include <cstddef>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "CharacterEventDiscrete.h"
#include "AdjacentRateModifier.h"
#include "CharacterHistoryRateModifier.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"

namespace RevBayesCore { class CharacterEvent; }

using namespace RevBayesCore;

AdjacentRateModifier::AdjacentRateModifier(size_t ns, size_t nc) : CharacterHistoryRateModifier(ns, nc),
    width(1),
    gain_factor(0.0),
    loss_factor(0.0),
    context_matrix( std::vector<std::vector<adjacency> >() ),
    context_set( std::vector<std::set<size_t> >() ),
    context_type("width"),
    exp_gain_factors( std::vector<double>(nc,1.0) ),
    exp_loss_factors( std::vector<double>(nc,1.0) )
{
    ;
}

AdjacentRateModifier::AdjacentRateModifier(const AdjacentRateModifier& g) : CharacterHistoryRateModifier(g)
{
    
    if (&g != this)
    {
        width = g.width;
        gain_factor = g.gain_factor;
        loss_factor = g.loss_factor;
        context_matrix = g.context_matrix;
        context_set = g.context_set;
        context_type = g.context_type;
        exp_gain_factors = g.exp_gain_factors;
        exp_loss_factors = g.exp_loss_factors;
    }
}

double AdjacentRateModifier::computeRateMultiplier(std::vector<CharacterEvent *> currState, CharacterEventDiscrete* newState, double age)
{
    
    
    if (context_type=="width")
    {
        return computeRateMultiplierUsingWidth(currState, newState, age);
    }
    else if (context_type=="matrix")
    {
        return computeRateMultiplierUsingMatrix(currState, newState, age);
    }
    else
    {
        std::stringstream ss;
        ss << "Undefined context_type (" << context_type << ")";
        throw RbException(ss.str());
    }
    
}


double AdjacentRateModifier::computeRateMultiplierUsingWidth(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age)
{
    std::vector<double> match(2, 0.0);
    
    size_t to_index = newState->getSiteIndex();
    size_t to_state = newState->getState();
    
    int lower_bound = ( (int)(to_index - width) < 0                     ? 0                               : (int)(to_index - width) );
    int upper_bound = ( (int)(to_index + width) >= this->num_characters ? (int)(this->num_characters - 1) : (int)(to_index + width) );
    
    for (int i = lower_bound; i <= upper_bound; i++)
    {
        if (i == to_index)
            continue;
        else if ( static_cast<CharacterEventDiscrete*>(currState[i])->getState() == to_state)
        {
            match[1] += 1.0;
        }
        else
        {
            match[0] += 1.0;
        }
        
    }
    
//    double r = exp(num_match * factor);
    double r = exp_gain_factors[ match[1] ] * exp_loss_factors[ match[0] ];
    //exp(gain_factor * match[1] + loss_factor * match[0]);
    
    return r;
}

double AdjacentRateModifier::computeRateMultiplierUsingMatrix(std::vector<CharacterEvent*> currState, CharacterEventDiscrete* newState, double age)
{
    std::vector<double> match(2, 0.0);
    
    size_t i = newState->getSiteIndex();
    size_t s = newState->getState();
    
    
    for (size_t j = 0; j < context_matrix[i].size(); j++)
    {
        const adjacency& edge = context_matrix[i][j];
        if (edge.from == edge.to)
        {
            continue;
        }
        else if ( static_cast<CharacterEventDiscrete*>(currState[edge.to])->getState() == s)
        {
            match[1] += edge.weight;
        }
        else
        {
            match[0] += edge.weight;
        }
    }
    
//    double r = exp(gain_factor * match[1] + loss_factor * match[0]);
    double r = exp_gain_factors[ match[1] ] * exp_loss_factors[ match[0] ];
    
    return r;
}


double AdjacentRateModifier::computeSiteRateMultiplier(const TopologyNode& node, CharacterEventDiscrete* currState, CharacterEventDiscrete* newState, double age)
{
    return 1.0;
}

double AdjacentRateModifier::computeSiteRateMultiplier(const TopologyNode& node, unsigned from, unsigned to, unsigned charIdx, double age)
{
    return 1.0;
}


AdjacentRateModifier* AdjacentRateModifier::clone(void) const
{
    return new AdjacentRateModifier(*this);
}

std::set<size_t> AdjacentRateModifier::getAffectedSites(CharacterEventDiscrete* newState) const
{
    std::set<size_t> s = context_set[ newState->getSiteIndex() ];
    
    return s;
}



void AdjacentRateModifier::setGainFactor(double f)
{
    // update factor
    gain_factor = f;
    
    // get exponent factor
    double m = exp(f);
    
    // reset factor vector
    exp_gain_factors = std::vector<double>( num_characters, 1.0 );
    exp_gain_factors[0] = 1.0;
    for (size_t i = 1; i < num_characters; i++)
    {
        exp_gain_factors[i] = exp_gain_factors[i-1] * m;
    }
}

void AdjacentRateModifier::setLossFactor(double f)
{
    // update factor
    loss_factor = f;
    
    // get exponent factor
    double m = exp(f);
    
    // reset factor vector
    exp_loss_factors = std::vector<double>( num_characters, 1.0 );
    exp_loss_factors[0] = 1.0;
    for (size_t i = 1; i < num_characters; i++)
    {
        exp_loss_factors[i] = exp_loss_factors[i-1] * m;
    }
}


void AdjacentRateModifier::setWidth(size_t w)
{
    context_type = "width";
    width = w;
}

void AdjacentRateModifier::setContextMatrix(const RbVector<RbVector<std::int64_t> >& c)
{

    context_type = "array";
    context_matrix = std::vector<std::vector<adjacency> >(this->num_characters);
    context_set = std::vector<std::set<size_t> >(this->num_characters);
    
    for (size_t i = 0; i < this->num_characters; i++)
    {
        for (size_t j = 0; j < this->num_characters; j++)
        {
            if (c[i][j] != 0.0)
            {
                // add adjacency
                adjacency v;
                v.from = i;
                v.to = j;
                v.weight = c[i][j];
                context_matrix[i].push_back(v);
                
                // add dependent context set -- site j whose rate depends on the state of site i
                context_set[i].insert(j);
            }
            
        }
        context_set[i].insert(i);
    }
}

void AdjacentRateModifier::setContextMatrix(const MatrixReal& c)
{
    
    context_type = "matrix";
    context_matrix = std::vector<std::vector<adjacency> >(this->num_characters);
    context_set = std::vector<std::set<size_t> >(this->num_characters);
    
    for (size_t i = 0; i < this->num_characters; i++)
    {
        for (size_t j = 0; j < this->num_characters; j++)
        {
            if (c[i][j] != 0.0)
            {
                // add adjacency
                adjacency v;
                v.from = i;
                v.to = j;
                v.weight = c[i][j];
                context_matrix[i].push_back(v);
                
                // add dependent context set -- site j whose rate depends on the state of site i
                context_set[i].insert(j);
            }
            
        }
        context_set[i].insert(i);
    }
}

void AdjacentRateModifier::update(void)
{
    ; // do nothing
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const AdjacentRateModifier& x)
{
    o << "AdjacentRateModifier";
    return o;
}
