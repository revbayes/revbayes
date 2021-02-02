#ifndef YAMLHelpRenderer_H
#define YAMLHelpRenderer_H

#include "RbHelpFunction.h"
#include "RbHelpType.h"
#include "YAML.h"

/**
 * @brief Class rendering the help to files in YAML format.
 *
 * The YAML-Help renderer converts the help entry into YAML files.
 *
 * \copyright (c) Copyright 2009-2014 (GPL version 3)
 * \author The RevBayes Development Core Team (Johan Dunfalk & Sebastian Hoehna)
 * \since Version 1.0, 2014-09-15
 */
namespace YAML
{
    class HelpRenderer
    {
        public:

            enum FunctionType { FUNCTION, METHOD, CONSTRUCTOR };

            std::string             renderHelp(const RevBayesCore::RbHelpType &typeHelp);
            std::string             renderHelp(const RevBayesCore::RbHelpFunction &functionHelp);

        private:

            std::vector<Element*> getAttributes(const RevBayesCore::RbHelpType &typeHelp);
            std::vector<Element*> getAttributes(const RevBayesCore::RbHelpFunction &functionHelp, FunctionType type = FUNCTION );
            std::vector<Element*> getAttributes(const RevBayesCore::RbHelpArgument &argHelp );
            std::vector<Element*> getAttributes(const RevBayesCore::RbHelpReference &refHelp);
    };
}
    

#endif

