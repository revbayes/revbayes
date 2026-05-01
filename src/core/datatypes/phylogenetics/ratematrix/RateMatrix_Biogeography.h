#ifndef RateMatrix_Biogeography_H
#define RateMatrix_Biogeography_H

#include "GeneralRateMatrix.h"
#include "TransitionProbabilityMatrix.h"
#include <complex>
#include <vector>
#include <map>
#include <list>


namespace RevBayesCore {
    
    class TransitionProbabilityMatrix;
    
    class RateMatrix_Biogeography : public GeneralRateMatrix {
        
    public:
        RateMatrix_Biogeography(size_t ns, size_t nc, size_t mrs);                                                                                               //!< Construct rate matrix with n states
        RateMatrix_Biogeography(const RateMatrix_Biogeography& m);                                                                                //!< Copy constructor
        virtual                         ~RateMatrix_Biogeography(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_Biogeography&                   operator=(const RateMatrix_Biogeography& r);
        
        // RateMatrix functions
        double                              averageRate(void) const;
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;   //!< Calculate the transition matrix
        RateMatrix_Biogeography*            clone(void) const;
        void                                fillRateMatrix(void);
        virtual RbVector<std::string>       getStateDescriptions(void) const;
        const RbVector<RbVector<double> >&  getDispersalRates(void) const;
        const RbVector<double>&             getExtirpationRates(void) const;                                                   //!< Return the extirpation rates
        virtual std::vector<double>         getStationaryFrequencies(void) const;                                              //!< Return the stationary frequencies

        void                                setDispersalRates(const RbVector<RbVector<double> >& dr);                          //!< Directly set dispersal rates
        void                                setExtirpationRates(const RbVector<double>& er);                        //!< Directly set extirpation rates
        
//        void                                setBirthRate(const double& br);
//        void                                setCladogeneticMatrix(const MatrixReal& cp);

        void                                update(void);
        
    private:
        std::string                         getRangeStr(const std::vector<unsigned>& v);
        void                                calculateCijk(void);                                                                //!< Do precalculations on
//        void                                computeConditionSurvival(MatrixReal& r);
//        void                                computeConditionSurvival(TransitionProbabilityMatrix& r);
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                initializeEigenSystem(void);                                                        //!< Create the eigen system, if it hasn't been created yet
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        
        size_t                          numBitsOn(std::vector<unsigned> v);
        void                            makeBits(void);
        void                            makeTransitions(void);
        std::string                     bitsToString( std::vector<unsigned> b ) const;
        std::string                     getStateDescriptions(void);
        
        std::vector<std::vector<unsigned> >                 bits;
        std::map<std::vector<unsigned>, unsigned>           inverseBits;
        std::vector<std::vector<std::vector<unsigned> > >   bitsByNumOn;
        std::vector<std::vector<unsigned> >                 statesToBitsByNumOn;
        std::map< std::vector<unsigned>, unsigned>          bitsToStatesByNumOn;
        std::vector<std::vector<unsigned> >                 transitions;
        std::vector<std::vector<unsigned> >                 lossOrGain;
        std::vector<std::vector<std::vector<unsigned> > >   transitionAreas;
        std::vector<std::vector<unsigned> >                 changedAreas;
        std::vector<std::vector<std::vector<unsigned> > >   affectingAreas;
        size_t                                              numCharacters;
        size_t                                              num_states;
        bool                                                useSquaring;
        
        bool                                                hasEigenSystem;                 //!< Whether we have initialized the eigen system yet
        EigenSystem*                                        theEigenSystem;                 //!< Holds the eigen system
        std::vector<double>                                 c_ijk;                          //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >                  cc_ijk;                         //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        // members
        RbVector<RbVector<double> >                         dispersalRates;
        RbVector<double>                                    extirpationRates;
        MatrixReal                                          cladogeneticMatrix;
        double                                              scalingFactor;
        double                                              birthRate;
        bool                                                useCladogenesis;
        bool                                                rescaleMatrix;
        size_t                                              maxRangeSize;
        TransitionProbabilityMatrix                         stationaryMatrix;
        
        mutable std::map<double,TransitionProbabilityMatrix>  storedTransitionProbabilities;
        mutable std::list<double>                             accessedTransitionProbabilities;
        unsigned                                              maxSizeStoredTransitionProbabilites;
        bool                                                  useStoredTransitionProbabilities;
    };
    
}

#endif /* defined(__revbayes_proj__RateMatrix_Biogeography__) */
