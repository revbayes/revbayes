/**
 * @file
 * This file contains the declaration of RnaState, which is
 * the class for the Rna data types in RevBayes.
 *
 * @brief Declaration of RnaState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-05-24 09:58:04 +0200 (Thu, 24 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RnaState.h 1568 2012-05-24 07:58:04Z hoehna $
 */

#ifndef RnaState_H
#define RnaState_H

#include <cstddef>
#include <ostream>

#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {

    class RnaState : public DiscreteCharacterState {
    
    public:
                                        RnaState(size_t n=4);                               //!< Default constructor
                                        RnaState(const std::string &s);                     //!< Constructor with nucleotide observation
        
        RnaState*                       clone(void) const;                                  //!< Get a copy of this object
    
        // Discrete character observation functions
        void                            addState(const std::string &symbol);                //!< Add a character state to the set of character states
        size_t                          getNumberOfStates(void) const;                      //!< Get the number of discrete states for the character
        RbBitSet                        getState(void) const;                               //!< Get the state (as the bitset)
        void                            setToFirstState(void);                              //!< Set this character state to the first (lowest) possible state
        void                            setStateByIndex(size_t index);                      //!< Set the discrete observation
        
        std::string                     getDataType(void) const;                            //!< Get the datatype as a common string.
        std::string                     getStateLabels(void) const;                         //!< Get valid state labels
        std::string                     getStringValue(void) const;                         //!< Get a representation of the character as a string
        void                            setState(const std::string &symbol);                //!< Compute the internal state value for this character.
        bool                            isGapState(void) const;                             //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const;                         //!< Get whether this is a missing character state
        void                            setGapState(bool tf);                               //!< set whether this is a gapped character
        void                            setMissingState(bool tf);                           //!< set whether this is a missing character
        
    private:
        
        char                            state;
    
    };
    
}

#endif

