#ifndef CharacterEventContinuous_H
#define CharacterEventContinuous_H

#include <cstdlib>
#include <vector>
#include <iosfwd>

#include "CharacterEvent.h"
#include "CharacterEventType.h"

namespace RevBayesCore {
    
    class CharacterEventContinuous : public CharacterEvent
    {
        
    public:
        CharacterEventContinuous(void);
        CharacterEventContinuous(size_t ch_idx, double val, double a, size_t t=CharacterEventType::UNDEFINED);
        CharacterEventContinuous(size_t ch_idx, const std::vector<double> &val, double a, size_t t=CharacterEventType::UNDEFINED);
        CharacterEventContinuous(const CharacterEventContinuous& c);
        ~CharacterEventContinuous(void);
        
        CharacterEventContinuous*           clone(void) const;
        double                              getState(size_t i) const;
        std::string                         getStateStr(void) const;
        void                                resize(size_t i);
        void                                setState(double s, size_t i);
        
        
    private:

        std::vector<double>                 value;                                             // the value/state of this event
        
        
    };
}


#endif

