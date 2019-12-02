#include "CharacterEventDiscrete.h"
#include "MatrixReal.h"
#include "RateGeneratorSequenceUsingMatrix.h"
#include "RateMatrix.h"
#include "RateMatrix_JC.h"
#include "RbVector.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateGeneratorSequenceUsingMatrix::RateGeneratorSequenceUsingMatrix(size_t ns, size_t nc) : RateGeneratorSequence(ns, nc)
{

    rateMatrix = new RateMatrix_JC(ns);
    rateModifiers = new RbVector<CharacterHistoryRateModifier>();

}



/** Copy constructor */
RateGeneratorSequenceUsingMatrix::RateGeneratorSequenceUsingMatrix(const RateGeneratorSequenceUsingMatrix& m) : RateGeneratorSequence( m )
{

    rateMatrix = m.rateMatrix->clone();
    rateModifiers = m.rateModifiers->clone();

}


/** Destructor */
RateGeneratorSequenceUsingMatrix::~RateGeneratorSequenceUsingMatrix(void)
{

    delete rateMatrix;
    delete rateModifiers;

}


RateGeneratorSequenceUsingMatrix& RateGeneratorSequenceUsingMatrix::operator=(const RateGeneratorSequenceUsingMatrix &r)
{

    RateGeneratorSequence::operator=( r );

    if (this != &r)
    {

        delete rateMatrix;
        delete rateModifiers;

        rateMatrix = r.rateMatrix->clone();
        rateModifiers = r.rateModifiers->clone();

    }

    return *this;
}



std::ostream& RevBayesCore::operator<<(std::ostream& o, const RateGeneratorSequenceUsingMatrix& x)
{
    std::streamsize previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();

    o.setf(previousFlags);
    o.precision(previousPrecision);

    return o;
}


//void RateGeneratorSequenceUsingMatrix::calculateTransitionProbabilities(TransitionProbabilityMatrix &P, double age) const
//{
//    const RateGenerator* rm = rateMatrix;
//
//    rm->calculateTransitionProbabilities(age, 0, 1.0, P);
//}



void RateGeneratorSequenceUsingMatrix::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    const RateGenerator* rm = rateMatrix;

    rm->calculateTransitionProbabilities(startAge, endAge, rate, P);
}


RateGeneratorSequenceUsingMatrix* RateGeneratorSequenceUsingMatrix::clone(void) const
{
    return new RateGeneratorSequenceUsingMatrix( *this );
}



std::set<size_t> RateGeneratorSequenceUsingMatrix::getAffectedSites(CharacterEventDiscrete* to) const
{
    size_t s_init[] = { to->getSiteIndex() };
    std::set<size_t> s( s_init, s_init+1 );

    for (size_t i = 0; i < rateModifiers->size(); i++) {
        std::set<size_t> tmp = (*rateModifiers)[i].getAffectedSites(to);
        s.insert( tmp.begin(), tmp.end() );
    }
    
    return s;
    
}

double RateGeneratorSequenceUsingMatrix::getRate(size_t from, size_t to, double age, double rate) const
{

    const RateGenerator* rm = rateMatrix;

    // no rate modifiers applied (per-site)
    double r = rm->getRate(from, to, age, rate);

    return r;

}


double RateGeneratorSequenceUsingMatrix::getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age, double rate) const
{
    size_t from_state = static_cast<CharacterEventDiscrete*>(from[ to->getSiteIndex() ])->getState();
    size_t to_state = to->getState();

    const RateGenerator* rm = rateMatrix;

    double r = rm->getRate(from_state, to_state, age, rate);
    
    for (size_t i = 0; i < rateModifiers->size(); i++)
    {
        r *= (*rateModifiers)[i].computeRateMultiplier(from, to, age);
    }

    return r;

}



double RateGeneratorSequenceUsingMatrix::getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, unsigned* counts, double age, double rate) const
{
    size_t from_state = static_cast<CharacterEventDiscrete*>(from[ to->getSiteIndex() ])->getState();
    size_t to_state = to->getState();
    
    const RateGenerator* rm = rateMatrix;
    
    double r = rm->getRate(from_state, to_state, age, rate);
    
    for (size_t i = 0; i < rateModifiers->size(); i++)
    {
        r *= (*rateModifiers)[i].computeRateMultiplier(from, to);
    }
    
    return r;

}

double RateGeneratorSequenceUsingMatrix::getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, std::vector<std::set<size_t> > sites_with_states, double age, double rate) const
{
    size_t from_state = static_cast<CharacterEventDiscrete*>(from[ to->getSiteIndex() ])->getState();
    size_t to_state = to->getState();
    
    const RateGenerator* rm = rateMatrix;
    
    double r = rm->getRate(from_state, to_state, age, rate);
    
    for (size_t i = 0; i < rateModifiers->size(); i++)
    {
        r *= (*rateModifiers)[i].computeRateMultiplier(from, to, sites_with_states, age);
    }
    
    return r;
}

double RateGeneratorSequenceUsingMatrix::getSiteRate(CharacterEventDiscrete* from, CharacterEventDiscrete* to, double age, double rate) const
{

    const RateGenerator* rm = rateMatrix;

    double r = rm->getRate(from->getState(), to->getState(), age, rate);
    
    return r;
}

double RateGeneratorSequenceUsingMatrix::getSiteRate(size_t from, size_t to, size_t charIdx, double age, double rate) const
{

    const RateGenerator* rm = rateMatrix;

    double r = rm->getRate(from, to, age, rate);

    throw(RbException("MJL: Will need to be fixed depending on how site rates are used."));
    
    return r;
}

std::vector<double> RateGeneratorSequenceUsingMatrix::getStationaryFrequencies(void) const
{
    const RateMatrix* rm = static_cast<RateMatrix*>(rateMatrix);
    if (rm != NULL) {
       return rm->getStationaryFrequencies();
    } else {
        // @MJL: Need to come up with a good way to manage stationary frequencies for e.g. biogeo process...
        throw(RbException("Model needs stationary frequencies!"));
    }

}

double RateGeneratorSequenceUsingMatrix::getSumOfRates(std::vector<CharacterEvent*> from, std::vector<size_t> counts, double age, double rate) const
{
    
    // get rate matrix
    const RateGenerator* rm = rateMatrix;
    
    // compute sum of leaving rates
    double sum = 0.0;
    
    // for each site in the starting sequence
    for (size_t i = 0; i < from.size(); i++)
    {
        size_t from_state = static_cast<CharacterEventDiscrete*>(from[i])->getState();
        
        // look at all outgoing states
        for (size_t to_state = 0; to_state < this->num_states; to_state++)
        {
            
            // ignore virtual events (where states match)
            if (from_state == to_state)
                continue;
            
            // get base rate
            double r = rm->getRate(from_state, to_state, age, 1.0);
            
            // get modified rate
            CharacterEventDiscrete to(i, to_state, age);
            for (size_t k = 0; k < rateModifiers->size(); k++)
            {
                CharacterHistoryRateModifier& chrm = (*rateModifiers)[k];
                double m = chrm.computeRateMultiplier(from, &to);
                r *= m;
            }
            
            // add rate to sum
            sum += r;
        }
    }
    
    sum *= rate;
    
    return sum;
}

double RateGeneratorSequenceUsingMatrix::getSumOfRates(std::vector<CharacterEvent*> from, double age, double rate) const
{

    std::vector<std::set<size_t> > sites_with_states(this->num_states);
    
    bool use_site_state_patterns = true;
    if (use_site_state_patterns)
    {
        for (size_t i = 0; i < from.size(); i++)
        {
            sites_with_states[ static_cast<CharacterEventDiscrete*>(from[i])->getState() ].insert(i);
        }
    }

    return getSumOfRates( from, sites_with_states, age, rate);
}

double RateGeneratorSequenceUsingMatrix::getSumOfRates(std::vector<CharacterEvent*> from, std::vector<std::set<size_t> > sites_with_states, double age, double rate) const
{
    
    // get rate matrix
    const RateGenerator* rm = rateMatrix;
    
    // compute sum of leaving rates
    double sum = 0.0;
    
    // for each site in the starting sequence
    for (size_t i = 0; i < from.size(); i++)
    {
        size_t from_state = static_cast<CharacterEventDiscrete*>(from[i])->getState();
        
        // look at all outgoing states
        for (size_t to_state = 0; to_state < this->num_states; to_state++)
        {
            
            // ignore virtual events (where states match)
            if (from_state == to_state)
                continue;
            
            // get base rate
            double r = rm->getRate(from_state, to_state, age, 1.0);
            
            // get modified rate
            CharacterEventDiscrete to(i, to_state, age);
            for (size_t k = 0; k < rateModifiers->size(); k++)
            {
                CharacterHistoryRateModifier& chrm = (*rateModifiers)[k];
                double m = chrm.computeRateMultiplier(from, &to, sites_with_states, age);
                r *= m;
                if (r == 0.0) break;
            }
            
            // add rate to sum
            sum += r;
        }
    }
    
    sum *= rate;
    
    return sum;

}


double RateGeneratorSequenceUsingMatrix::getSumOfRatesDifferential(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, std::vector<std::set<size_t> > sites_with_states, double age, double rate) const
{
    double r = 0.0;
    
    size_t index = to->getSiteIndex();
    size_t old_state = static_cast<CharacterEventDiscrete*>(from[ index ])->getState();
    size_t new_state = to->getState();
    
    CharacterEventDiscrete* possible_event = new CharacterEventDiscrete(*to);
    
    std::set<size_t> affected_sites = getAffectedSites(to);
    
//    std::vector<std::set<size_t> > old_sites_with_states = sites_with_states;
    std::vector<std::set<size_t> > new_sites_with_states = sites_with_states;
    new_sites_with_states[ old_state ].erase( index );
    new_sites_with_states[ new_state ].insert( index );
    
    for (size_t s = 0; s < num_states; s++)
    {
        possible_event->setState(s);
        
        // subtract the contribution of rates leaving the old state
        if (s != old_state)
        {
            r -= getRate(from, possible_event, sites_with_states, age, rate);
        }
        
        // prepare for new state
        static_cast<CharacterEventDiscrete*>(from[ index ])->setState( new_state );
        
        // add the contribution of rates leaving the new state
        if (s != new_state)
        {
            r += getRate(from, possible_event, new_sites_with_states, age, rate);
        }
        
        // revert to old state
        static_cast<CharacterEventDiscrete*>(from[ index ])->setState( old_state );
    }
    
    delete possible_event;
    
    return r;
}

double RateGeneratorSequenceUsingMatrix::getSumOfRatesDifferential(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age, double rate) const
{
    // this will return the rate differential
    double r = 0.0;

    
    // get the old and new states for the transition
    size_t to_index = to->getSiteIndex();
    size_t old_to_state = static_cast<CharacterEventDiscrete*>(from[ to->getSiteIndex() ])->getState();
    size_t new_to_state = to->getState();
    
    // collect the sites affected by the transition
    std::set<size_t> affected_sites = getAffectedSites(to);
    
    // subtract the sum of rates contribution leaving the old from state
    for (std::set<size_t>::iterator it = affected_sites.begin(); it != affected_sites.end(); it++) {
        
        // index of affected site
        size_t index = *it;
        
        // create a dummy event for the affected site
        CharacterEventDiscrete* possible_event = new CharacterEventDiscrete( *(static_cast<CharacterEventDiscrete*>(from[index])) );
        
        // first subtract the sum of rates leaving each old site-state
        size_t old_from_state = static_cast<CharacterEventDiscrete*>(from[index])->getState();
        for (size_t s = 0; s < num_states; s++) {
            if (s != old_from_state)
            {
                possible_event->setState(s);
                r -= getRate(from, possible_event, age, rate);
            }
        }
        
        // apply the new to-state
        static_cast<CharacterEventDiscrete*>(from[to_index])->setState(new_to_state);
        
        // then add the sum of rates leaving each new site-state
        size_t new_from_state = static_cast<CharacterEventDiscrete*>(from[index])->getState();
        for (size_t s = 0; s < num_states; s++) {
            if (s != new_from_state)
            {
                possible_event->setState(s);
                r += getRate(from, possible_event, age, rate);
            }
        }
        
        // restore the old to-state
        static_cast<CharacterEventDiscrete*>(from[to_index])->setState(old_to_state);
        
        // free memory for the dummy state
        delete possible_event;
    }
    
    // return the rate differential
    return r;
}

void RateGeneratorSequenceUsingMatrix::setRateMatrix(const RateGenerator* r)
{
    if (r != rateMatrix)
    {
        delete rateMatrix;
        rateMatrix = r->clone();
    }
}


void RateGeneratorSequenceUsingMatrix::setRateModifiers(const RbVector<CharacterHistoryRateModifier>* r)
{
    if (r != rateModifiers)
    {
        delete rateModifiers;
        rateModifiers = r->clone();
    }
}


void RateGeneratorSequenceUsingMatrix::updateMap(void)
{
    if (needsUpdate)
    {
        ; // do nothing ...

        needsUpdate = false;
    }
}
