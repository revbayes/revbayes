/**
 * @file
 * This file contains the declaration of DnaState, which is
 * the class for the DNA data types in RevBayes.
 *
 * @brief Declaration of DnaState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-05-24 09:58:04 +0200 (Thu, 24 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: DnaState.h 1568 2012-05-24 07:58:04Z hoehna $
 */

#ifndef DnaState_H
#define DnaState_H

#include <cstddef>
#include <ostream>

#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {

    class DnaState : public DiscreteCharacterState {
    
    public:
                                        DnaState(size_t n=4);                               //!< Default constructor
                                        DnaState(const std::string &s);                     //!< Constructor with nucleotide observation
                                        DnaState(const RbBitSet& bs);                       //!< Constructor with which letters are observed.

        bool                            operator==(const CharacterState& x) const;
        bool                            operator==(const DnaState& x) const;
//        bool                            operator<(const DnaState& x) const;
        void                            operator+=(int);
        void                            operator-=(int);

        bool                            isAmbiguous(void) const;
        size_t                          getStateIndex(void) const;                          //!< Get the index of the current state

        DnaState*                       clone(void) const override;                         //!< Get a copy of this object

        // Discrete character observation functions

        void                            addState(const std::string &symbol) override;       //!< Add a character state to the set of character states
        size_t                          getNumberOfStates(void) const override;             //!< Get the number of discrete states for the character
        RbBitSet                        getState(void) const override;                      //!< Get the state (as the bitset)
        void                            setToFirstState(void) override;                     //!< Set this character state to the first (lowest) possible state
        void                            setStateByIndex(size_t index) override;             //!< Set the discrete observation
        std::string                     getDataType(void) const override;                   //!< Get the datatype as a common string.
        std::string                     getStateLabels(void) const override;                //!< Get valid state labels
        std::string                     getStringValue(void) const override;                //!< Get a representation of the character as a string
        void                            setState(const std::string &symbol) override;       //!< Compute the internal state value for this character.
        bool                            isGapState(void) const override;                    //!< Get whether this is a gapped character state
        bool                            isMissingState(void) const override;                //!< Get whether this is a missing character state
        void                            setGapState(bool tf) override;                      //!< set whether this is a gapped character
        void                            setMissingState(bool tf) override;                  //!< set whether this is a missing character
        
    private:
        
        char                            state;
    
    };
    
}

#endif

