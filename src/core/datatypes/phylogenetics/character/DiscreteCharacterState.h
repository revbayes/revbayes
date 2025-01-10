#ifndef DiscreteCharacterState_H
#define DiscreteCharacterState_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "CharacterState.h"
#include "RbBitSet.h"

namespace RevBayesCore {


    /** @brief Abstract base class for discrete characters data types.
     *
     * To accommodate ambiguous characters, derived classes can implement the current state
     * as a set of states instead of a single state.
    **/
    class DiscreteCharacterState : public CharacterState {

    public:
        virtual                                ~DiscreteCharacterState(void) {}

        virtual bool                            operator==(const CharacterState& x) const;
        bool                                    operator!=(const CharacterState& x) const;
        virtual bool                            operator<(const CharacterState& x) const;
        void                                    operator++();       //!< Increment
        void                                    operator++(int i);  //!< Increment
        virtual void                            operator+=(int i);  //!< Increment
        void                                    operator--();       //!< Decrement
        void                                    operator--(int i);  //!< Decrement
        virtual void                            operator-=(int i);  //!< Decrement

        virtual bool                            isAmbiguous(void) const;
        virtual std::string                     getStringValue(void) const;
        virtual std::string                     getStateDescription(void) const;        //!< Get a description of the current state
        virtual std::vector<std::string>        getStateDescriptions(void) const;       //!< Get all possible state labels as a vector
        virtual size_t                          getNumberObservedStates(void) const;    //!< How many states are in the set of current states
        virtual size_t                          getNumberOfStates(void) const;
        virtual size_t                          getStateIndex(void) const;              //!< Get the index of the current state
        virtual bool                            isStateSet(size_t index) const;         //!< Is this state part of the current set?

        virtual DiscreteCharacterState*         clone(void) const = 0;

        virtual void                            addState(const std::string &symbol) = 0;    //!< Add a state to the set of current states
        virtual RbBitSet                        getState(void) const = 0;                   //!< Get the current state (as a bitset)
        virtual void                            setToFirstState(void) = 0;                  //!< Set the state to the first (lowest) possible state
        virtual void                            setStateByIndex(size_t index) = 0;          //!< Set the current state (by index)
        virtual void                            setState(const std::string &symbol) = 0;    //!< Set the current state (by state)

        // Discrete character observation functions
        virtual std::string                     getStateLabels(void) const = 0;             //!< Get labels for all possible states
        const std::vector<double>&              getWeights() const;                         //!< Get the weights of each state
        bool                                    isWeighted() const;                         //!< Is the state weighted?
        void                                    setWeighted( bool tf );                     //!< Sets whether the state is weighted or not

    protected:
                                                DiscreteCharacterState(size_t n);   //!< Constructor
    };

}

#endif
