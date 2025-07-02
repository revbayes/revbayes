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
        virtual                                            ~MarginalLikelihoodEstimator(void);                                         //!< Virtual destructor
        
        // public methods
        virtual MarginalLikelihoodEstimator*                clone(void) const = 0;                                                     //!< Create a new deep copy
        virtual double                                      marginalLikelihood( void ) const = 0;
        double                                              getESS(const std::vector<double> values) const;                            //!< Calculate the effective sample size for an arbitrary vector of numeric values
        
        // block bootstrap functions
        size_t                                              offsetModulo(size_t i, size_t n);                                          //!< Helper function
        std::vector<size_t>                                 getIndices(std::pair<size_t, size_t> a, size_t n);                         //!< Get a vector of indices to select given a start index and vector length
        
    protected:
        
        // members
        std::vector< double >                               powers;
        std::vector< std::vector< double> >                 likelihoodSamples;
    };
    
}

#endif
