#ifndef RateGeneratorSequence_Biogeography_H
#define RateGeneratorSequence_Biogeography_H

#include "RateGeneratorSequenceUsingMatrix.h"
#include "DistanceRateModifier.h"
#include <vector>

namespace RevBayesCore {
    
    class BranchHistory;
    class RateGeneratorSequence_Biogeography : public RateGeneratorSequenceUsingMatrix {
        
    public:
        RateGeneratorSequence_Biogeography(size_t nc, bool fe=true, unsigned mrs=0);                                                                               //!< Construct rate matrix with n states
        RateGeneratorSequence_Biogeography(const RateGeneratorSequence_Biogeography& m);                                                           //!< Copy constructor
        virtual                                     ~RateGeneratorSequence_Biogeography(void);                                       //!< Destructor
        
        // overloaded operators
        RateGeneratorSequence_Biogeography&                       operator=(const RateGeneratorSequence_Biogeography& r);
        
        // virtual RateGeneratorSequence functions
//        void                                        calculateTransitionProbabilities(TransitionProbabilityMatrix& P, double age=0.0) const;   //!< Calculate the transition probabilities for the rate matrix
        void                                        calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrixmatrix
//        void                                        calculateTransitionProbabilities(TransitionProbabilityMatrix& P, double startAge, double endAge, double rate, size_t charIdx) const;   //!< Calculate the transition matrixmatrix
        RateGeneratorSequence_Biogeography*                       clone(void) const;
        double                                      getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, double age=0.0, double rate=1.0) const;
        double                                      getRate(std::vector<CharacterEvent*> from, CharacterEventDiscrete* to, unsigned* counts, double age=0.0, double rate=1.0) const;
        double                                      getRate(size_t from, size_t to, double age=0.0, double rate=1.0) const;
        double                                      getSiteRate(CharacterEventDiscrete* from, CharacterEventDiscrete* to, double age=0.0, double rate=1.0) const;
        double                                      getSiteRate(size_t from, size_t to, size_t charIdx=0, double age=0.0, double rate=1.0) const;
        double                                      getSumOfRates(std::vector<CharacterEvent*> from, double age=0.0, double rate=1.0) const;
        double                                      getSumOfRates(std::vector<CharacterEvent*> from, std::vector<size_t> counts, double age=0.0, double rate=1.0) const;
        double                                      getUnnormalizedSumOfRates(std::vector<CharacterEvent*> from, std::vector<size_t> counts, double age=0.0, double rate=1.0) const;
        const bool                                  isAreaAvailable(size_t charIdx, double age=0.0) const;
        const bool                                  areAreasAdjacent(size_t fromCharIdx, size_t toCharIdx, double age=0.0) const;
        const std::set<size_t>                      getRangeAndFrontierSet(BranchHistory* bh, double age=0.0) const;
        void                                        updateMap(void);
        
        // public methods
        double                                      getDistancePower(void) const;
        void                                        setDistancePower(double d);
        
        // other crazy stuff for BiogeographyRateGeneratorSequenceFunction to handle
        void                                        setDistanceRateModifier(const DistanceRateModifier& gdrm);
        void                                        setGeographicDistancePowers(const DistanceRateModifier& gdrm);
        const DistanceRateModifier &                getDistanceRateModifier(void);
        const std::vector<double>&                  getEpochs(void) const;
        
    private:
        size_t                                      numOn(const std::vector<CharacterEvent*>& s) const;
        size_t                                      numOn(const std::vector<CharacterEvent*>& s, double age) const;
        unsigned                                    getEpochIndex(double age) const;
        size_t                                      branchOffset;
        double                                      distancePower;
        unsigned                                    maxRangeSize;
        
        // geography models
        bool                                        useDistanceRateModifier;
        DistanceRateModifier*                       distanceRateModifier;
        std::vector<double>                         epochs;
        size_t                                      numEpochs;
        size_t                                      epochOffset;
        std::vector<double>                         adjacentAreaVector;
        std::vector<double>                         availableAreaVector;

        // model flags
        bool                                        forbidExtinction;
        
        bool                                        useUnnormalizedRates;
        bool                                        useAreaAdjacency;
        bool                                        useAreaAvailable;
        bool                                        useDistanceDependence;
        bool                                        useRootFrequencies;
        
    };
    
}
#endif
