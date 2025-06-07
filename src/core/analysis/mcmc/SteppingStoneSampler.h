#ifndef SteppingStoneSampler_H
#define SteppingStoneSampler_H

#include <iosfwd>

#include "MarginalLikelihoodEstimator.h"

namespace RevBayesCore {
    
    /**
     * @brief SteppingStone-Sampler class.
     *
     * The SteppingStone sampler analyzes the output of a power posterior analysis
     * and computes the SteppingStone-Sampler marginal likelihood.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2012-06-17
     *
     */
    class SteppingStoneSampler : public MarginalLikelihoodEstimator {
        
    public:
        SteppingStoneSampler(const std::string &fn, const std::string &pn, const std::string &ln, const std::string &del);             //!< Constructor initializing the object.
        virtual                                            ~SteppingStoneSampler(void);                                                //!< Virtual destructor
        
        // public methods
        SteppingStoneSampler*                               clone(void) const;                                                         //!< Create a deep copy
        double                                              marginalLikelihood( void ) const;                                          //!< Compute the marginal likelihood using stepping-stone sampling
        double                                              getESS(const std::vector<double> values) const;                            //!< Calculate the effective sample size for an arbitrary vector of numeric values
        double                                              standardError( void ) const;                                               //!< Compute the standard error of the marginal likelihood estimate
    };
    
}

#endif
