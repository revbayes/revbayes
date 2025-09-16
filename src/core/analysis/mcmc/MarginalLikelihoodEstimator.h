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
     * @author The RevBayes Development Core Team (Sebastian Hoehna, David Cerny)
     * @since Version 1.0, 2014-06-24
     *
     */
    class MarginalLikelihoodEstimator : public Cloneable, public Parallelizable {
        
    public:
        MarginalLikelihoodEstimator(const path &fn, const std::string &pn, const std::string &ln, const std::string &del);
        virtual                                            ~MarginalLikelihoodEstimator(void);                                         //!< Virtual destructor
        
        // public methods
        virtual MarginalLikelihoodEstimator*                clone(void) const = 0;                                                     //!< Create a new deep copy
        virtual double                                      marginalLikelihoodGeneral(std::vector<double> power_vec, std::vector< std::vector<double> > lnl_vec) const = 0;
        double                                              marginalLikelihood( void ) const;                                          //!< Compute the marginal likelihood for this analysis
        virtual double                                      standardError( void ) const = 0;                                           //!< Compute the standard error for this analysis
        
        // block bootstrap functions
        std::int64_t                                        offsetModulo(std::int64_t i, std::int64_t n) const;                        //!< Helper function
        std::vector<std::int64_t>                           getIndices(std::pair<std::int64_t, std::int64_t> a, std::int64_t n) const; //!< Get a vector of indices to select given a start index and vector length
        std::vector< std::vector< std::vector<double> > >   blockBootstrap(size_t repnum, double prop, bool print_to_file) const;      //!< Generate block bootstrap replicates
        double                                              standardErrorBlockBootstrap(size_t repnum, double prop, bool print, bool verbose) const; //!< Compute the standard error from block bootstrap replicates
        
    protected:
        
        // members
        path                                                filename;
        std::vector<double>                                 powers;
        std::vector< std::vector<double> >                  likelihoodSamples;
    };
    
}

#endif
