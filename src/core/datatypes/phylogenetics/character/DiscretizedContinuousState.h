/**
 * @file
 * This file contains the declaration of DiscretizedContinuousState, which is
 * the base class for character data types that are represented
 * as discretized versions of continuous data in RevBayes.
 *
 * @brief Declaration of DiscretizedContinuousState
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 */


#ifndef DiscretizedContinuousState_H
#define DiscretizedContinuousState_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {
    
    class DiscretizedContinuousState : public DiscreteCharacterState {
        
    public:
        DiscretizedContinuousState(size_t n=10, double d=0.01);                                         //!< Default constructor
        DiscretizedContinuousState(std::vector<std::string> labels, std::vector<double> pts, double d); //!< Constructor with state labels
        DiscretizedContinuousState(const std::string &s, std::vector<double> pts, int m, double d);     //!< Constructor with an observation
        DiscretizedContinuousState(int s, std::vector<double> pts, int m, double d);                    //!< Constructor with an observation
                
        DiscretizedContinuousState*            clone(void) const;                                  //!< Get a copy of this object
        
        // Discrete character observation functions
        void                            addState(const std::string &symbol);                //!< Add a character state to the set of character states
        RbBitSet                        getState(void) const;                               //!< Get the state (as the bitset)
        void                            setToFirstState(void);                              //!< Set this character state to the first (lowest) possible state
        void                            setState(const std::string &symbol);                //!< Compute the internal state value for this character.
        void                            setStateByIndex(size_t index);                      //!< Set the discrete observation

        void                            addState(int s);                                    //!< Add the state with the given index.
        void                            addStateDescriptions(const std::vector<std::string>& d);
        std::string                     getDataType(void) const;                            //!< Get the datatype as a common string.
        double                          getDeltaX(void) const;
        std::vector<double>             getPoints(void) const;
        std::string                     getStateDescription(void) const;
        std::vector<std::string>        getStateDescriptions(void) const;
        std::string                     getStateLabels(void) const;                         //!< Get valid state labels
        std::string                     getStringValue(void) const;                         //!< Get a representation of the character as a string
        bool                            isGapState(void) const;                             //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const;                         //!< Get whether this is a missing character state
        void                            setGapState(bool tf);                               //!< set whether this is a gapped character
        void                            setMissingState(bool tf);                           //!< set whether this is a missing character
        
        const std::vector<double>&      getWeights( void ) const;                            //!< Get the weight of the state
        bool                            isWeighted( void ) const;
        void                            setWeighted( bool tf );
        void                            setWeights(std::vector<double> w);                  //!< set the weight for each character
        
    private:
        
        bool                            is_gap = false;
        bool                            is_missing = false;
        size_t                          index_single_state = 0;
        size_t                          num_observed_states = 0;
        RbBitSet                        state;
        std::vector<std::string>        state_descriptions;
        double                          dx;
        std::vector<double>             points;
        std::vector<double>             weights;                                            //!< Weights are used when the "average" option is used
        bool                            weighted;


    };

}

#endif
