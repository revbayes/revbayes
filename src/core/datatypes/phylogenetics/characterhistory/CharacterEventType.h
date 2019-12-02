//
//  CharacterEventType.h
//  revbayes-proj
//
//  Created by Michael Landis on 2/6/17.
//  Copyright © 2017 Michael Landis. All rights reserved.
//

#ifndef CharacterEventType_h
#define CharacterEventType_h

namespace RevBayesCore {
    namespace CharacterEventType {
        const size_t NUM_EVENTS          = 3;
        
        const size_t UNDEFINED           = 0;
        const size_t STATE_CHANGE        = 1;
        const size_t CLASS_CHANGE        = 2;
    };
};



#endif /* CharacterEventType_h */
