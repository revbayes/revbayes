#ifndef RateGeneratorSequenceUsingMatrix_H
#define RateGeneratorSequenceUsingMatrix_H

#include "CharacterHistoryRateModifier.h"
#include "RateGeneratorSequence.h"
#include "RateMatrix.h"
#include "RbVector.h"
#include "TransitionProbabilityMatrix.h"

#include <set>
#include <vector>

namespace RevBayesCore {
    
    class RateGeneratorSequenceUsingMatrix : public RateGeneratorSequence {
        
    public:
        RateGeneratorSequenceUsingMatrix(const RateGeneratorSequenceUsingMatrix& m);                                                              //!< Copy constructor
        RateGeneratorSequenceUsingMatrix(size_t ns, size_t nc);                                                          //!< Construct rate matrix with n states
        RateGeneratorSequenceUsingMatrix&                 operator=(const RateGeneratorSequenceUsingMatrix& r);
        
        virtual                            ~RateGeneratorSequenceUsingMatrix(void);                                                                 //!< Destructor
        
        
        // virtual public methods
//        virtual void                        calculateTransitionProbabilities(TransitionProbabilityMatrix& P, double age=0.0) const;   //!< Calculate the transition probabilities for the rate matrix
        virtual void                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrixmatrix
//        virtual void                        calculateTransitionProbabilities(TransitionProbabilityMatrix& P, size_t charIdx) const;
        virtual RateGeneratorSequenceUsingMatrix*         clone(void) const;
        
        virtual std::set<size_t>            getAffectedSites(CharacterEventDiscrete* to) const;
        
        virtual double                      getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, unsigned* counts, double age=0.0, double rate=1.0) const;
        virtual double                      getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age=0.0, double rate=1.0) const;
       
        virtual double                      getRate(size_t from, size_t to, double age=0.0, double rate=1.0) const;
        virtual double                      getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, std::vector<std::set<size_t> > sites_with_states, double age=0.0, double rate=1.0) const;
        
        virtual double                      getSiteRate(CharacterEventDiscrete* from, CharacterEventDiscrete* to, double age=0.0, double rate=1.0) const;
        virtual double                      getSiteRate(size_t from, size_t to, size_t charIdx=0, double age=0.0, double rate=1.0) const;
        
        virtual std::vector<double>         getStationaryFrequencies(void) const;
        
        virtual double                      getSumOfRates(std::vector<CharacterEvent*> from, double age=0.0, double rate=1.0) const;
        virtual double                      getSumOfRates(std::vector<CharacterEvent*> from, std::vector<size_t> counts, double age=0.0, double rate=1.0) const;
        virtual double                      getSumOfRates(std::vector<CharacterEvent*> from, std::vector<std::set<size_t> > sites_with_states, double age=0.0, double rate=1.0) const;
        
        virtual double                      getSumOfRatesDifferential(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age=0.0, double rate=1.0) const;
        virtual double                      getSumOfRatesDifferential(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to,  std::vector<std::set<size_t> > sites_with_states, double age=0.0, double rate=1.0) const;
        
        void                                setRateMatrix(const RateGenerator *rm);
        void                                setRateModifiers(const RbVector<CharacterHistoryRateModifier>* rm);
        virtual void                        updateMap(void);
        
        
    protected:
        // protected members available for derived classes
        RateGenerator*                          rateMatrix;
        RbVector<CharacterHistoryRateModifier>* rateModifiers;
        
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const RateGeneratorSequenceUsingMatrix& x);                                           //!< Overloaded output operator
    
}

#endif
