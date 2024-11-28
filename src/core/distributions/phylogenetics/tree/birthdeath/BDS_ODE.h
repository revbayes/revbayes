#ifndef BDS_ODE_H
#define BDS_ODE_H

#include "AbstractBirthDeathProcess.h"
#include "RateMatrix.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Birth-Death-Shift differential stepper equation
     *
     */
    class BDS_ODE {
        
    public:
        
        BDS_ODE( const std::vector<double> &m, const RateGenerator* q, double r, bool backward_time, bool extinction_only, bool allow_shifts_extint=true );
        
        void operator() ( const std::vector< double > &x, std::vector< double > &dxdt , const double t );
        
        void            setEventMap( const std::map<std::vector<unsigned>, double> &e );
        void            setSpeciationRate( const std::vector<double> &s );
        void            setSerialSamplingRate( const std::vector<double> &s );
        
    private:
        
        std::vector<double>                         mu;                                 //!< vector of extinction rates, one rate for each character state
        std::vector<double>                         lambda;                             //!< vector of speciation rates, one rate for each character state
        std::vector<double>                         psi;                                //!< vector of fossilization rates, one rate for each character state
        size_t                                      num_states;                         //!< the number of character states = q->getNumberOfStates()
        const RateGenerator*                        Q;                                  //!< anagenetic rate matrix
        std::map<std::vector<unsigned>, double>     event_map;                          //!< cladogenetic event map, with the structure pair< [ancestor_state, daughter_1_state, daughter_2_state], speciation_rate >
        double                                      rate;                               //!< clock rate for anagenetic change
        
        // flags to modify behabior
        bool                                        extinction_only;                    //!< calculate only extinction probabilities
        bool                                        use_speciation_from_event_map;      //!< do we use the speciation rates from the event map?
        bool                                        backward_time;                      //!< computation backward in time (otherwise forward)?
        bool                                        allow_rate_shifts_extinction;       //!< should we allow for rate-shifts in the extinct lineages?

    };
    
}

#endif
