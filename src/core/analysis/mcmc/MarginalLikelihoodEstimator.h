#ifndef MarginalLikelihoodEstimator_H
#define MarginalLikelihoodEstimator_H


#include <vector>
#include <iosfwd>

#include "Cloneable.h"
#include "Parallelizable.h"
#include "RbFileManager.h"

namespace RevBayesCore {
    
    /**
     * @brief Marginal likelihood estimator interface.
     *
     * This interface provides some common methods needed for a path-sampler
     * and a stepping-stone-sampler.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2014-06-24
     *
     */
    class MarginalLikelihoodEstimator : public Cloneable, public Parallelizable {
        
    public:
        MarginalLikelihoodEstimator(const path &fn, const std::string &pn, const std::string &ln, const std::string &del);
        virtual                                            ~MarginalLikelihoodEstimator(void);                                                          //!< Virtual destructor
        
        // public methods
        virtual MarginalLikelihoodEstimator*                clone(void) const = 0;                                                                      //!< Create a new deep copy
        virtual double                                      marginalLikelihood( void ) const = 0;
        
    protected:
        
        // members
        std::vector< double >                               powers;
        std::vector< std::vector< double> >                 likelihoodSamples;
    };
    
}

#endif
