#ifndef MultiOptionRule_H
#define MultiOptionRule_H

#include <string>
#include <vector>

#include "OptionRule.h"

namespace RevLanguage {
class Argument;
class RevObject;
class RlString;

/**
 * An OptionRule that allows several options at once, given as a vector of strings or as one
 * string with the options joined by a separator, e.g. "noabsencesites|nosingletonpresence"
 * (the MrBayes style). Each option may appear at most once. The argument always reaches the
 * function as a String[].
 */
class MultiOptionRule : public OptionRule {

    public:
        MultiOptionRule(const std::string& argName, RlString* defVal, const std::vector<std::string>& optVals, const std::string& sep, const std::string& argDesc );

        MultiOptionRule*                    clone(void) const;
        Argument                            fitArgument(Argument& arg) const;
        double                              isArgumentValid(Argument &arg) const;

    private:
        bool                                selected(const RevObject& value, std::vector<std::string>& tokens) const;   //!< The options in value, false if it is not a string or string vector

        std::string                         separator;
};

}

#endif
