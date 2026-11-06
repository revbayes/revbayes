/**
 * @file
 * This file contains the declaration of Delimiter, which is
 * used to describe a dot dot dot formal (...), representing
 * a variable number of arguments. Delimiter must always be
 * used as the last argument rule. It will match any number of
 * arguments, including 0. Types are not checked.
 *
 * @brief Declaration of Delimiter
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 *
 * $Id$
 */

#ifndef Delimiter_H
#define Delimiter_H

#include "ArgumentRule.h"

namespace RevLanguage {


class Delimiter : public ArgumentRule {

    static constexpr const char* DESCRIPTION = "The field separator string. Values on each line of the file are separated by this string. If sep = \"\" the separator is \'white space\', that is one or more spaces, tabs, newlines or carriage returns.";

public:
    Delimiter( const std::string& def = "", const std::string &desc = DESCRIPTION);         //! Some type specification needs to be met, value arguments by default

    // Basic utility functions
    Delimiter* clone(void) const { return new Delimiter(*this); }   //!< Clone object
};


// When writing a delimiter, we cannot use patterns like "any white" space that we use to split fields.
// Instead we must use a specific string.
class WriteDelimiter : public ArgumentRule {

    static constexpr const char* DESCRIPTION = "The field separator string. It may not be empty.";

public:
    WriteDelimiter( const std::string& def, const std::string &desc = DESCRIPTION);        //! Some type specification needs to be met, value arguments by default

    // Basic utility functions
    WriteDelimiter* clone(void) const { return new WriteDelimiter(*this); }   //!< Clone object
};

}

#endif

