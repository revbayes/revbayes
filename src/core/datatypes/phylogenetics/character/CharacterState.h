#ifndef Character_H
#define Character_H

#include <cstddef>
#include <ostream>

#include "Cloneable.h"

namespace RevBayesCore {

    /** @brief Abstract base class for character data types (RNA, DNA, etc.) in RevBayes
    **/
    class CharacterState : public Cloneable {
    
    public:
        virtual                        ~CharacterState(void) {}
        
        virtual bool                    operator==(const CharacterState& x) const = 0;      //!< Equality
        virtual bool                    operator!=(const CharacterState& x) const = 0;      //!< Inequality
        virtual bool                    operator<(const CharacterState& x) const = 0;       //!< Less than
    
        // pure virtual
        virtual CharacterState*         clone(void) const = 0;
        virtual std::string             getDataType(void) const = 0;                        //!< Get the datatype as a common string, e.g. DNA, RNA or Protein
        virtual size_t                  getNumberOfStates(void) const = 0;                  //!< Get the number of states for this character
        virtual std::string             getStringValue(void) const = 0;                     //!< Get a representation of the character as string
        virtual bool                    isAmbiguous(void) const = 0;                        //!< Is the character missing or ambiguous
        virtual bool                    isGapState(void) const = 0;                         //!< Get whether this is a gapped character state
        virtual bool                    isMissingState(void) const = 0;                     //!< Get whether this is a missing character state
        virtual void                    setGapState(bool tf) = 0;                           //!< Set whether this is a gapped character state
        virtual void                    setMissingState(bool tf) = 0;                       //!< Set whether this is a missing character state

    protected:
                                        CharacterState(void);
        
    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const CharacterState& x);  //!< Overloaded output operator

}

#endif

