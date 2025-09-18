/**
 * @file
 * This file contains the declaration of Expanded Standard State, which *is the class for the expanded standard data type in RevBayes.
 *
 * @brief Declaration of ExpandedStandardState
 *
 * (c) Copyright 2025-
 * @author The RevBayes Development Team
 *
 */
#ifndef ExpandedStandardState_H
#define ExpandedStandardState_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "DiscreteCharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore
{

  class ExpandedStandardState : public DiscreteCharacterState
  {

  public:
    ExpandedStandardState(size_t n = 10);
    ExpandedStandardState(const std::string &s, int m);
    ExpandedStandardState(int s, int m);

    ExpandedStandardState *clone(void) const; //!< Get a copy of this object

    // Discrete character observation functions
    void addState(const std::string &symbol); //!< Add a character state to the set of character states
    RbBitSet getState(void) const;            //!< Get the state (as the bitset)
    void setToFirstState(void);               //!< Set this character state to the first (lowest) possible state
    void setState(const std::string &symbol); //!< Compute the internal state value for this character.
    void setStateByIndex(size_t index);       //!< Set the discrete observation

    void addState(int s); //!< Add the state with the given index.
    void addStateDescriptions(const std::vector<std::string> &d);
    std::string getDataType(void) const; //!< Get the datatype as a common string.
    std::string getStateDescription(void) const;
    std::vector<std::string> getStateDescriptions(void) const;
    std::string getStateLabels(void) const;    //!< Get valid state labels
    std::string getStringValue(void) const;    //!< Get a representation of the character as a string
    bool isGapState(void) const;               //!< Get whether this is a gapped character state
    bool isMissingState(void) const;           //!< Get whether this is a missing character state
    void setGapState(bool tf);                 //!< set whether this is a gapped character
    void setMissingState(bool tf);             //!< set whether this is a missing character
    void setStateLabels(const std::string &l); //!< set the state labels

  private:
    bool is_gap;
    bool is_missing;
    size_t index_single_state;
    size_t num_observed_states;
    RbBitSet state;
    std::vector<std::string> state_descriptions;
    std::string labels;
  };
}

#endif
