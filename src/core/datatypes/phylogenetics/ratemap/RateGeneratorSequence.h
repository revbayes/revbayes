#ifndef RateGeneratorSequence_H
#define RateGeneratorSequence_H

#include "CharacterEvent.h"
#include "CharacterEventDiscrete.h"
#include "CharacterState.h"
#include "Cloneable.h"
#include "RateGenerator.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"

#include <set>
#include <vector>

namespace RevBayesCore {
    
//    class RateGeneratorSequence : public RateGenerator, public Cloneable {
    class RateGeneratorSequence : public RateGenerator {
    
    public:
        RateGeneratorSequence(size_t ns, size_t nc);                                                          //!< Construct rate matrix with n states
        
        virtual                            ~RateGeneratorSequence(void);                                                                 //!< Destructor
        
        // public methods
        size_t                              getNumberOfCharacters(void) const;                  //!< Return the number of characters

        // virtual public methods
        virtual void                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const = 0;   //!< Calculate the transition matrixmatrix
        virtual RateGeneratorSequence*      clone(void) const = 0;
        
        virtual double                      getRate(size_t from, size_t to, double age=0.0, double rate=1.0) const = 0;
        virtual double                      getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age=0.0, double rate=1.0) const = 0;
//        virtual double                      getRate(std::vector<CharacterEvent*> from, CharacterEvent* to, unsigned* counts, double rate=1.0, double age=0.0) const = 0;
        virtual double                      getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, std::vector<std::set<size_t> > sites_with_states, double age=0.0, double rate=1.0) const = 0;
        
        virtual std::vector<double>         getStationaryFrequencies(void) const = 0;
        
        using RateGenerator::getSumOfRates;
        virtual double                      getSumOfRates(std::vector<CharacterEvent*> from, std::vector<std::set<size_t> > sites_with_states, double age=0.0, double rate=1.0) const = 0;
        
        virtual double                      getSumOfRatesDifferential(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age=0.0, double rate=1.0) const = 0;
        virtual double                      getSumOfRatesDifferential(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, std::vector<std::set<size_t> > sites_with_states, double age=0.0, double rate=1.0) const = 0;
      

        virtual void                        updateMap(void);
        
    protected:
        // protected members available for derived classes
//        size_t                              numStates;                                          //!< The number of character states
        size_t                              num_characters;                                      //!< The number of characters
        bool                                needsUpdate;
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const RateGeneratorSequence& x);                                           //!< Overloaded output operator
    
}

#endif
