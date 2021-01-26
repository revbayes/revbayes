/**
 * @file
 * This file contains the declaration of NaturalNumbersState, which is
 * the base class for character data types that are represented
 * as natural numbers (such as morphological characters) in RevBayes.
 * Gap states represented as -1
 *
 * @brief Declaration of NaturalNumbersState
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 */


#ifndef NaturalNumbersState_H
#define NaturalNumbersState_H

#include <stddef.h>
#include <ostream>
#include <set>
#include <vector>

#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {
    
    class NaturalNumbersState : public DiscreteCharacterState {
        
    public:
        NaturalNumbersState(size_t s=0, size_t m=0);                                        //!< Constructor with an observation
        NaturalNumbersState(const std::string &s, size_t m=0);                              //!< Constructor with an observation
                
        NaturalNumbersState*            clone(void) const;                                  //!< Get a copy of this object
        
        // Discrete character observation functions
        void                            addState(const std::string &symbol);                //!< Add a character state to the set of character states
        RbBitSet                        getState(void) const;                               //!< Get the state (as the bitset)
        void                            setToFirstState(void);                              //!< Set this character state to the first (lowest) possible state
        void                            setState(const std::string &symbol);                //!< Compute the internal state value for this character.
        void                            setStateByIndex(size_t index);                      //!< Set the discrete observation

        void                            addState(size_t s);                                 //!< Add the state with the given index.
        void                            addStateDescriptions(const std::vector<std::string>& d);
        std::string                     getDataType(void) const;                            //!< Get the datatype as a common string.
        std::string                     getStateDescription(void) const;
        std::vector<std::string>        getStateDescriptions(void) const;
        std::string                     getStateLabels(void) const;                         //!< Get valid state labels
        std::string                     getStringValue(void) const;                         //!< Get a representation of the character as a string
        bool                            isGapState(void) const;                             //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const;                         //!< Get whether this is a missing character state
        bool                            isPositiveState(void) const;                        //!< Get whether this is a positive character state
        void                            setGapState(bool tf);                               //!< set whether this is a gapped character
        void                            setMissingState(bool tf);                           //!< set whether this is a missing character
        void                            setPositiveState(bool tf);                          //!< set whether this is a positive character
        
    private:
        
        bool                            is_gap;
        bool                            is_missing;
        bool                            is_positive;
        size_t                          single_state;
        size_t                          num_observed_states;
        size_t                          max_state;
        std::set<size_t>                state;
        std::vector<std::string>        state_descriptions;

    };
}

#endif
