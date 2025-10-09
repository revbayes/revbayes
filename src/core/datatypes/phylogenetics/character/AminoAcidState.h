#ifndef AminoAcidState_H
#define AminoAcidState_H

#include <cstddef>
#include <ostream>

#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {

    
    
    /**
     * Amino Acid State class.
     *
     * Amino acid states are represented by their own data structure (this one).
     * Internally we store the amino acid state as an unsigned integer number.
     * We use the first 20 bits of this variable to mark which state is set.
     * This allows us to use ambiguous observations.
     * We provide the common functions, see DiscreteCharacterState.h
     *
     * \copyright (c) Copyright 2009- (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-06-21
     */
    class AminoAcidState : public DiscreteCharacterState {
    
    public:
                                        AminoAcidState(size_t n=20);                        //!< Default constructor
                                        AminoAcidState(const std::string &s);               //!< Constructor with nucleotide observation
        
        AminoAcidState*                 clone(void) const;                                  //!< Get a copy of this object

        // Discrete character observation functions
        void                            addState(const std::string &symbol);                //!< Add a character state to the set of character states
        RbBitSet                        getState(void) const;                               //!< Get the state (as the bitset)
        void                            setToFirstState(void);                              //!< Set this character state to the first (lowest) possible state
        void                            setState(const std::string &symbol);                //!< Compute the internal state value for this character.
        void                            setStateByIndex(size_t index);                      //!< Set the discrete observation

        std::string                     getDataType(void) const;                            //!< Get the datatype as a common string.
        std::string                     getStateLabels(void) const;                         //!< Get valid state labels
        bool                            isGapState(void) const;                             //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const;                         //!< Get whether this is a missing character state
        void                            setGapState(bool tf);                               //!< set whether this is a gapped character
        void                            setMissingState(bool tf);                           //!< set whether this is a missing character

    private:
        
        bool                            is_gap;
        bool                            is_missing;
        size_t                          index_single_state;
        size_t                          num_observed_states;
        RbBitSet                        state;

    };
    
}

#endif

