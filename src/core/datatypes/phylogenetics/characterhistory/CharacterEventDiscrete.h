#ifndef CharacterEventDiscrete_H
#define CharacterEventDiscrete_H

#include <cstdlib>
#include <iosfwd>

#include "CharacterEvent.h"
#include "CharacterEventType.h"

namespace RevBayesCore {
    
    class CharacterEventDiscrete : public CharacterEvent
    {

    public:
        CharacterEventDiscrete(void);
        CharacterEventDiscrete(size_t i, size_t s, double a, size_t t=CharacterEventType::UNDEFINED);
        CharacterEventDiscrete(const CharacterEventDiscrete& c);
        ~CharacterEventDiscrete(void);

        CharacterEventDiscrete*             clone(void) const;
        size_t                              getState(void) const;
        std::string                         getStateStr(void) const;
        void                                setState(size_t s);

    protected:


    private:

        size_t                              state;                                             // from statespace -- possibly make a set of state/positions...
        

    };
}


#endif
