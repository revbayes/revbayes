#ifndef OdeHeterogeneousRateBirthDeath_H
#define OdeHeterogeneousRateBirthDeath_H

#include <cstddef>
#include <vector>

#include "RbVector.h"

//typedef boost::array< double , 4 > state_type;
//typedef std::vector< double > state_type;

namespace RevBayesCore {
    
    /**
     * @brief Multi-rate birth-death ODE.
     *
     *
     * This class represents the ordinary differential equations for the
     * multi-rate birth-death process.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-03-11, version 1.0
     *
     */
    class OdeHeterogeneousRateBirthDeath {
        
    public:
        
        OdeHeterogeneousRateBirthDeath( const RbVector<double> &l, const RbVector<double> &m, double r, bool a);
        
        void operator() ( const std::vector< double > &x , std::vector< double > &dxdt , const double t );

        void                        setCurrentRateCategory(size_t i);
        
    private:
        
        RbVector<double>            lambda;
        RbVector<double>            mu;
        double                      switch_rate;
        size_t                      num_categories;
        size_t                      current_rate_category;
        bool                        allow_same_category;
        
    };
    
}

#endif
