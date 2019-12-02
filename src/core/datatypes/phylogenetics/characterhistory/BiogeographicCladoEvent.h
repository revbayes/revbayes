#ifndef BiogeographicCladoEvent_H
#define BiogeographicCladoEvent_H

namespace RevBayesCore {
    namespace BiogeographicCladoEvent {
        const unsigned NUM_STATES          = 6;
        
        const unsigned SYMPATRY_NARROW     = 0;         // A    -> A    | A
        const unsigned SYMPATRY_WIDESPREAD = 1;         // ABCD -> ABCD | ABCD
        const unsigned SYMPATRY_SUBSET     = 2;         // ABCD -> ABCD | A
        const unsigned ALLOPATRY           = 3;         // ABCD -> AB   |   CD
        const unsigned PARAPATRY           = 4;         // ABCD -> ABC  |  BCD
        const unsigned JUMP_DISPERSAL      = 5;         // ABC  -> ABC  |    D
    };
};


#endif /* defined(BiogeographicCladoEvent_H) */
