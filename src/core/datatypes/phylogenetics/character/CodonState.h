/**
 * @file
 * This file contains the declaration of CodonState, which is
 * the class for the Codon data types in RevBayes.
 *
 * @brief Declaration of CodonState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-05-24 09:58:04 +0200 (Thu, 24 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: CodonState.h 1568 2012-05-24 07:58:04Z hoehna $
 */

#ifndef CodonState_H
#define CodonState_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "AminoAcidState.h"
#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {
    
    class CodonState : public DiscreteCharacterState {
        
    public:
        CodonState(size_t n=0);                                                             //!< Default constructor
        CodonState(const std::string &symbol);                                              //!< Constructor with codon observation
        
        CodonState*                     clone(void) const;                                  //!< Get a copy of this object
        
        // Discrete character observation functions
        void                            addState(const std::string &symbol);                //!< Add a character state to the set of character states
        RbBitSet                        getState(void) const;                               //!< Get the state (as the bitset)
        void                            setToFirstState(void);                              //!< Set this character state to the first (lowest) possible state
        void                            setState(const std::string &symbol);                //!< Compute the internal state value for this character.
        void                            setStateByIndex(size_t index);                      //!< Set the discrete observation

        AminoAcidState                  getAminoAcidState() const;
        std::string                     getDataType(void) const;                            //!< Get the datatype as a common string.
        std::string                     getStateLabels(void) const;                         //!< Get valid state labels
        std::string                     getStringValue(void) const;                         //!< Get a representation of the character as a string
        std::vector<unsigned int>       getTripletStates() const;
        bool                            isStopCodon(void) const;                            //!< Is the codon a stop codon
        bool                            isGapState(void) const;                             //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const;                         //!< Get whether this is a missing character state
        void                            setGapState(bool tf);                               //!< set whether this is a gapped character
        void                            setMissingState(bool tf);                           //!< set whether this is a missing character
        
        static std::string              CODONS[];
        static size_t                   CODON_TO_TRIPLET_INDICES[];
        static size_t                   TRIPLET_TO_CODON_INDICES[];

    private:
        
        bool                            is_gap;
        bool                            is_missing;
        size_t                          index_single_state;
        size_t                          num_observed_states;
        RbBitSet                        state;

        
    };
    
}

#endif

