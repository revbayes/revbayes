#ifndef CDCladoSE_H
#define CDCladoSE_H

#include "AbstractBirthDeathProcess.h"
#include "RateMatrix.h"

#include <vector>
#include <boost/array.hpp>

typedef std::vector< double > state_type;

namespace RevBayesCore {
    
    /**
     * @brief Multi-rate cladogenetic birth-death ODE.
     *
     *
     * This class represents the ordinary differential equations for the character dependent
     * cladogenetic multi-rate birth-death process (ClaSSE: Goldberg and Igic, 2012)
     * Will Freyman 6/22/16
     *
     */
    class CDCladoSE {
        
    public:
        
        CDCladoSE( const std::vector<double> &m, const RateGenerator* q, std::map<std::vector<unsigned>, double> e, double r);
        
        void operator() ( const state_type &x , state_type &dxdt , const double t );
        
    private:
        
        std::vector<double>                         mu;             // vector of extinction rates, one rate for each character state
        size_t                                      num_states;     // the number of character states = q->getNumberOfStates()
        const RateGenerator*                        Q;              // anagenetic rate matrix
        std::map<std::vector<unsigned>, double>     eventMap;       // cladogenetic event map, with the structure pair< [ancestor_state, daughter_1_state, daughter_2_state], speciation_rate >
        double                                      rate;           // clock rate for anagenetic change
        
        
    };
    
}

#endif
