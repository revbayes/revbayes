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
     * @author The RevBayes Development Core Team (Sebastian Hoehna, David Cerny)
     * @since Version 1.0, 2012-06-17
     *
     */
    class SteppingStoneSampler : public MarginalLikelihoodEstimator {
        
    public:
        SteppingStoneSampler(const std::string &fn, const std::string &pn, const std::string &ln, const std::string &del);             //!< Constructor initializing the object.
        virtual                                            ~SteppingStoneSampler(void);                                                //!< Virtual destructor
        
        // public methods
        SteppingStoneSampler*                               clone(void) const;                                                         //!< Create a deep copy
        double                                              marginalLikelihoodGeneral(std::vector<double> power_vec, std::vector< std::vector<double> > lnl_vec) const; //!< Compute the stepping-stone marginal likelihood for arbitrary values
        double                                              standardError(void) const;                                                 //!< Compute the standard error of the marginal likelihood estimate
    };
    
}

#endif
