#ifndef DiscreteCharacterState_H
#define DiscreteCharacterState_H

#include "CharacterState.h"
#include "RbBitSet.h"

#include <ostream>
#include <string>
#include <vector>

namespace RevBayesCore {

    class DiscreteCharacterState : public CharacterState {

    public:
        virtual                                ~DiscreteCharacterState(void) {}

        bool                                    operator==(const CharacterState& x) const;          //!< Equality
        bool                                    operator!=(const CharacterState& x) const;          //!< Inequality
        bool                                    operator<(const CharacterState& x) const;           //!< Less than
        void                                    operator++();                                       //!< Increment
        void                                    operator++(int i);                                  //!< Increment
        void                                    operator+=(int i);                                  //!< Increment
        void                                    operator--();                                       //!< Decrement
        void                                    operator--(int i);                                  //!< Decrement
        void                                    operator-=(int i);                                  //!< Decrement

        bool                                    isAmbiguous(void) const;                            //!< Is the character missing or ambiguous
        virtual std::string                     getStringValue(void) const;                         //!< Get a representation of the character as a string
        virtual std::string                     getStateDescription(void) const;
        virtual std::vector<std::string>        getStateDescriptions(void) const;
        virtual size_t                          getNumberObservedStates(void) const;                //!< How many states are observed for the character
        virtual size_t                          getNumberOfStates(void) const;                      //!< Get the number of discrete states for the character
        virtual size_t                          getStateIndex(void) const;                          //!< Get the state index
        virtual bool                            isStateSet(size_t index) const;                     //!< Is the discrete observation


        virtual DiscreteCharacterState*         clone(void) const = 0;                              //!< Get a copy of this object


        virtual void                            addState(const std::string &symbol) = 0;                //!< Add a character state to the set of character states
        virtual RbBitSet                        getState(void) const = 0;                               //!< Get the state (as the bitset)
        virtual void                            setToFirstState(void) = 0;                              //!< Set this character state to the first (lowest) possible state
        virtual void                            setStateByIndex(size_t index) = 0;                      //!< Set the discrete observation
        virtual void                            setState(const std::string &symbol) = 0;                //!< Set the discrete observation

        // Discrete character observation functions
        virtual std::string                     getStateLabels(void) const = 0;                     //!< Get valid state labels
        const std::vector<double>&              getWeights() const;                                //!< Get the weights of the state
        bool                                    isWeighted() const;                                //!< Is the state weighted?
        void                                    setWeighted( bool tf );                                //!< Sets whether the state is weighted or not

    protected:
                                                DiscreteCharacterState(size_t n);                   //!< Constructor
        
//        size_t                                  num_states;
        bool                                    weighted;
        std::vector<double>                     weights;
    };

}

#endif
