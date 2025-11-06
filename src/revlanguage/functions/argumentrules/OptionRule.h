/**
 * @file
 * This file contains the declaration of OptionRule, which is
 * used to describe argument rules corresponding to the
 * selection of one of several RlString options.
 *
 * @brief Declaration of OptionRule
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef OptionRule_H
#define OptionRule_H

#include <iosfwd>
#include <vector>

#include "ArgumentRule.h"

namespace RevLanguage {
class Argument;
class RlString;

class OptionRule : public ArgumentRule {

    public:
        OptionRule(const std::string& argName, const std::vector<std::string>& optVals, const std::string& argDesc );                                        //!< Constructor of rule without default value
        OptionRule(const std::string& argName, RlString* defVal, const std::vector<std::string>& optVals, const std::string& argDesc );                      //!< Constructor of rule with default value

        // Basic utility functions
        OptionRule*                         clone(void) const;                                                                  //!< Clone object
        const std::vector<std::string>&     getOptions(void) const;                                                             //!< Get the options
        virtual double                      isArgumentValid(Argument &arg) const;                                               //!< Is var a valid argument?
        void                                printValue(std::ostream& o) const;                                                  //!< Print value for user

    protected:
        bool                                areOptionsUnique(const std::vector<std::string>& optVals) const;                    //!< Test if options are unique

        std::vector<std::string>            options;                                                                            //!< Permissible values
    
};

}

#endif

