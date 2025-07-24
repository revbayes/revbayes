/**
 * @file
 * This file contains the declaration of the RevBayes main file which runs the RevLanguage environment.
 *
 *
 * @brief Declaration of the main
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @since Version 1.0, 2012-08-02
 *
 * $Id$
 */


#ifndef RevLanguageMain_H
#define RevLanguageMain_H

#include <string>
#include <vector>
#include <optional>

class RevLanguageMain {
    
public:
    
    RevLanguageMain(bool interactive, bool echo, bool error_exit, bool quiet);
    
    void startRevLanguageEnvironment(const std::vector<std::string> &expressions, const std::optional<std::string>& filename, const std::vector<std::string> &args);

private:
    
    bool interactive;           // keep asking stdin for commands until q()
    bool echo;                  // print commands to the screen
    bool error_exit;            // quit on the first error
    bool quiet;                 // suppress header
};

#endif
