#ifndef ConditionalPosteriorOrdinate_H
#define ConditionalPosteriorOrdinate_H

#include <iosfwd>
#include <vector>

#include "Cloneable.h"
#include "RbFileManager.h"
#include "Parallelizable.h"


namespace RevBayesCore {
    
    /**
     * @brief ConditionalPosteriorOrdinate (CPO) class.
     *
     * The CPO analyzes the output of an MCMC run to compute the leave-one-out cross-validation fitness.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.2, 2022-09-12
     *
     */
    class ConditionalPosteriorOrdinate : public Cloneable, public Parallelizable {
        
    public:
        ConditionalPosteriorOrdinate(const path& fn, const std::string& del, const std::vector<std::string>& skip_col_names);           //!< Constructor initializing the object.
        virtual                                            ~ConditionalPosteriorOrdinate(void);                                         //!< Virtual destructor
        
        // public methods
        ConditionalPosteriorOrdinate*                       clone(void) const;                                                          //!< Create a deep copy
        double                                              predictiveProbability( const std::vector<double>& counts, bool site_probs_as_log, bool folded ) const;                                           //!< Compute the marginal likelihood using SteppingStone-Sampler
        
    private:
        
        std::vector< std::vector<double> >                  samples;
    };
    
}

#endif
