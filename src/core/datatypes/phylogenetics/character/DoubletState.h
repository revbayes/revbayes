/**
 * @file
 * This file contains the declaration of DoubletState, which is
 * the class for the Doublet data types in RevBayes.
 *
 * @brief Declaration of DoubletState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-05-24 09:58:04 +0200 (Thu, 24 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: DoubletState.h 1568 2012-05-24 07:58:04Z hoehna $
 */

#ifndef DoubletState_H
#define DoubletState_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {

    class DoubletState : public DiscreteCharacterState {

    public:
        DoubletState(size_t n=0);                                                             //!< Default constructor
        DoubletState(const std::string &symbol);                                              //!< Constructor with codon observation

        DoubletState*                     clone(void) const;                                  //!< Get a copy of this object

        // Discrete character observation functions
        void                            addState(const std::string &symbol);                //!< Add a character state to the set of character states
        RbBitSet                        getState(void) const;                               //!< Get the state (as the bitset)
        void                            setToFirstState(void);                              //!< Set this character state to the first (lowest) possible state
        void                            setState(const std::string &symbol);                //!< Compute the internal state value for this character.
        void                            setStateByIndex(size_t index);                      //!< Set the discrete observation

        std::string                     getDataType(void) const;                            //!< Get the datatype as a common string.
        std::string                     getStateLabels(void) const;                         //!< Get valid state labels
        std::string                     getStringValue(void) const;                         //!< Get a representation of the character as a string
        std::vector<unsigned int>       getDoubletStates() const;

        bool                            isGapState(void) const;                             //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const;                         //!< Get whether this is a missing character state
        void                            setGapState(bool tf);                               //!< set whether this is a gapped character
        void                            setMissingState(bool tf);                           //!< set whether this is a missing character

        static std::string              DOUBLETS[];

    private:

        bool                            is_gap;
        bool                            is_missing;
        size_t                          index_single_state;
        size_t                          num_observed_states;
        RbBitSet                        state;
    };
    
}

#endif

