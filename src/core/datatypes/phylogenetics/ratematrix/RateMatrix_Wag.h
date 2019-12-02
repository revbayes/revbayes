///**
// * @file
// * This file contains the declaration of RateMatrix_JC, which is
// * class that holds a rate matrix for a continuous-time Markov model.
// *
// * @brief Declaration of RateMatrix_JC
// *
// * (c) Copyright 2009-
// * @date Last modified: $Date: 2012-12-05 12:20:37 -0800 (Wed, 05 Dec 2012) $
// * @author The RevBayes Development Core Team
// * @license GPL version 3
// *
// * $Id: RateMatrix.h 1901 2012-12-05 20:20:37Z boussau $
// */
//
//#ifndef RateMatrix_Wag_H
//#define RateMatrix_Wag_H
//
//#include "RateMatrix_Empirical.h"
//
//
//namespace RevBayesCore {
//    
//    class TransitionProbabilityMatrix;
//    
//    class RateMatrix_Wag : public RateMatrix_Empirical {
//        
//    public:
//        RateMatrix_Wag(void);                                                                                                //!< Construct rate matrix with n states
//        virtual                             ~RateMatrix_Wag(void);                                                               //!< Destructor
//        
//        // RateMatrix functions
//        RateMatrix_Wag*                     clone(void) const;
//        
//    private:
//        
//        
//    };
//    
//}
//
//#endif
//

#ifndef RateMatrix_Wag_H
#define RateMatrix_Wag_H

#include "TimeReversibleRateMatrix.h"
#include "EigenSystem.h"
#include <complex>
#include <vector>


namespace RevBayesCore {
    
    //    class EigenSystem;
    class TransitionProbabilityMatrix;
    
    /**
     * @brief WAG (Whelan and Goldman 2001) rate matrix class.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-07-04, version 1.0
     */
    class RateMatrix_Wag : public TimeReversibleRateMatrix {
        
    public:
        RateMatrix_Wag(void);                                                                                                   //!< Construct rate matrix with n states
        RateMatrix_Wag(const RateMatrix_Wag& m);                                                                                //!< Copy constructor
        virtual                             ~RateMatrix_Wag(void);                                                              //!< Destructor
        
        // overloaded operators
        RateMatrix_Wag&                     operator=(const RateMatrix_Wag& r);
        
        // RateMatrix functions
        virtual RateMatrix_Wag&             assign(const Assignable &m);                                                                                            //!< Assign operation that can be called on a base class instance.
        void                                calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const;    //!< Calculate the transition matrix
        RateMatrix_Wag*                     clone(void) const;
        void                                update(void);
        virtual void                        initFromString( const std::string &s );                                             //!< Serialize (resurrect) the object from a string value
        
    private:
        void                                calculateCijk(void);                                                                //!< Do precalculations on eigenvectors and their inverse
        void                                setEmpiricalStationaryFrequencies(void);
        void                                setEmpiricalExchangeabilityRates(void);
        void                                tiProbsEigens(double t, TransitionProbabilityMatrix& P) const;                      //!< Calculate transition probabilities for real case
        void                                tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const;               //!< Calculate transition probabilities for complex case
        void                                updateEigenSystem(void);                                                            //!< Update the system of eigenvalues and eigenvectors
        
        EigenSystem*                        theEigenSystem;                                                                     //!< Holds the eigen system
        std::vector<double>                 c_ijk;                                                                              //!< Vector of precalculated product of eigenvectors and their inverse
        std::vector<std::complex<double> >  cc_ijk;                                                                             //!< Vector of precalculated product of eigenvectors and thier inverse for complex case
        
        
    };
    
}

#endif


