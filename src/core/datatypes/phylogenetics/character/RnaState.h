/**
 * @file
 * This file contains the declaration of RnaState, which is
 * the class for the RNA data types in RevBayes.
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
                                        RnaState(const RbBitSet& bs);                       //!< Constructor with which letters are observed.

        bool                            operator==(const CharacterState& x) const override;
        bool                            operator==(const RnaState& x) const;
        void                            operator+=(int) override;
        void                            operator-=(int) override;
        
        RnaState*                       clone(void) const override;                         //!< Get a copy of this object
    
        // Discrete character observation functions
        void                            addState(const std::string &symbol) override;       //!< Add a character state to the set of character states
        std::string                     getDataType(void) const override;                   //!< Get the datatype as a common string.
        size_t                          getNumberOfStates(void) const override;             //!< Get the number of discrete states for the character
        RbBitSet                        getState(void) const override;                      //!< Get the state (as the bitset)
        size_t                          getStateIndex(void) const override;                 //!< Get the index of the current state
        std::string                     getStateLabels(void) const override;                //!< Get valid state labels
        std::string                     getStringValue(void) const override;                //!< Get a representation of the character as a string
        bool                            isAmbiguous(void) const override;
        bool                            isGapState(void) const override;                    //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const override;                //!< Get whether this is a missing character state
        void                            setGapState(bool tf) override;                      //!< set whether this is a gapped character
        void                            setMissingState(bool tf) override;                  //!< set whether this is a missing character
        void                            setState(const std::string &symbol) override;       //!< Compute the internal state value for this character.
        void                            setStateByIndex(size_t index) override;             //!< Set the discrete observation
        void                            setToFirstState(void) override;                     //!< Set this character state to the first (lowest) possible state
        
    private:
        
        char                            state;
    
    };

    inline bool operator==(const DiscreteCharacterState& lhs, const RnaState& rhs)
    {
        return rhs.operator==(lhs);
    }

    inline bool operator==(const RnaState& rhs, const DiscreteCharacterState& lhs)
    {
        return lhs.operator==(rhs);
    }
    
}

#endif

