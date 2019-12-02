#ifndef MultispeciesCoalescentMigrationODE_H
#define MultispeciesCoalescentMigrationODE_H

#include "AbstractBirthDeathProcess.h"
#include "RateMatrix.h"

#include <vector>
#include <boost/array.hpp>

typedef std::vector< double > state_type;

namespace RevBayesCore {
    
    /**
     * @brief Multispecies-coalescent-migration ODE.
     *
     *
     * This class represents the ordinary differential equations for the
     * multispecies-coalescent-migration model.
     * Sebastian Hoehna 8/11/17
     *
     */
    class MultispeciesCoalescentMigrationODE {
        
        public:
        
        MultispeciesCoalescentMigrationODE( const std::vector<double> &t, const RateGenerator* q, double r, size_t n_ind, size_t n_pop);
        
        void operator() ( const state_type &x , state_type &dxdt , const double t );
                
        private:
        
        std::vector<double>                         theta;                              //!< vector of population sizes
        size_t                                      num_individuals;                    //!< the number of individuals
        size_t                                      num_populations;                    //!< the number of individuals
        const RateGenerator*                        Q;                                  //!< migration rate matrix
        double                                      rate;                               //!< clock rate for migration rates
        
        
    };
    
}

#endif
