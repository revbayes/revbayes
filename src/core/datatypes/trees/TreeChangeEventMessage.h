//
//  TreeChangeEventMessage.hpp
//  revbayes-proj
//
//  Created by Michael Landis on 12/4/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef TreeChangeEventMessage_hpp
#define TreeChangeEventMessage_hpp

namespace RevBayesCore {
    namespace TreeChangeEventMessage {
        const unsigned DEFAULT           = 0;
        const unsigned TOPOLOGY          = 1; // e.g. resulting from NNI move
        const unsigned BRANCH_LENGTH     = 2; // e.g. resulting from branch scale update
        const unsigned CHARACTER_HISTORY = 3; // e.g. resulting from mvCharacterHistory
    };
};

#endif /* TreeChangeEventMessage_hpp */
