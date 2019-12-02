//
//  MpiUtilities.h
//  revbayes-proj
//
//  Created by Michael Landis on 4/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__MpiUtilities__
#define __revbayes_proj__MpiUtilities__

#include <sstream>
#include <string>

namespace RevBayesCore {
    
    namespace MpiUtilities {
        
        void DebugWait(int rank);
        void DebugMsg(const std::stringstream& s);
        void DebugMsg(const std::string& s);
        void DebugMsg(const std::string& s, int x);
        void DebugMsg(const std::string& s, double x);
        void DebugMsgPid(const std::string& s, int pid);


    }
}

#endif /* defined(__revbayes_proj__MpiUtilities__) */

