/**
 * @file
 * This file contains the declaration of Character, which is
 * the abstract base class for characters, such as DNA, RNA, 
 * etc. data types in RevBayes.
 *
 * @brief Declaration of Character
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-05-24 09:58:04 +0200 (Thu, 24 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: Character.h 1568 2012-05-24 07:58:04Z hoehna $
 */

#ifndef Character_H
#define Character_H

#include "Cloneable.h"

#include <ostream>
#include <string>

namespace RevBayesCore {

    class CharacterState : public Cloneable {
    
    public:
        virtual                        ~CharacterState(void) {}
        
        virtual bool                    operator==(const CharacterState& x) const = 0;      //!< Equality
        virtual bool                    operator!=(const CharacterState& x) const = 0;      //!< Inequality
        virtual bool                    operator<(const CharacterState& x) const = 0;       //!< Less than
    
        // pure virtual
        virtual CharacterState*         clone(void) const = 0;                              //!< Get a copy of this object
        virtual std::string             getDataType(void) const = 0;                        //!< Get the datatype as a common string, e.g. DNA, RNA or Protein
        virtual size_t                  getNumberOfStates(void) const = 0;                  //!< Get the number states for a character (return 0 in base class)
        virtual std::string             getStringValue(void) const = 0;                     //!< Get a representation of the character as string
        virtual bool                    isAmbiguous(void) const = 0;                        //!< Is the character missing or ambiguous
        virtual bool                    isGapState(void) const = 0;                         //!< Get whether this is a gapped character state
        virtual bool                    isMissingState(void) const = 0;                     //!< Get whether this is a missing character state
        virtual void                    setGapState(bool tf) = 0;                           //!< set whether this is a gapped character
        virtual void                    setMissingState(bool tf) = 0;                       //!< set whether this is a missing character

    protected:
                                        CharacterState(void);                               //!< Default constructor
        
    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const CharacterState& x);                                //!< Overloaded output operator

}

#endif

