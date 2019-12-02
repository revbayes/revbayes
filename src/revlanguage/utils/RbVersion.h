/**
 * @file
 * This file contains the declaration of RbVersion,
 * which holds and prints information about the current RevBayes version
 *
 * @brief Declaration of RbVersion
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since version 1.0 2012-02-25
 *
 * $Id$
 */

#ifndef RbVersion_H
#define RbVersion_H

#include <string>

class RbVersion {
    
public:
    RbVersion(void);                                                            //!< Default constructor
    
    // Access functions
    std::string                     getDate(void) const;
    std::string                     getHeader(void) const;
    std::string                     getGitBranch(void) const;
    std::string                     getGitCommit(void) const;
    std::string                     getVersion(void) const;                     //!< Retrieve the version
    void                            printHeader(std::ostream& o);
    
    
};

#endif
