#ifndef CharacterEvent_H
#define CharacterEvent_H

#include <cstdlib>
#include <iosfwd>

#include "Cloneable.h"
#include "CharacterEventType.h"

namespace RevBayesCore {
    
    class CharacterEvent : public Cloneable
    {
        
    public:
        ~CharacterEvent(void);
        
        // pure virtual methods
        virtual CharacterEvent*             clone(void) const = 0;
        virtual std::string                 getStateStr(void) const = 0;

        
        // public methods implemented for derived classes
        
        // getters
        double                              getAge(void) const;
        size_t                              getEventType(void) const;
        size_t                              getSiteIndex(void) const;
        
        // setters
        void                                setAge(double a);
        void                                setEventType(size_t t);
        void                                setMissingState(bool tf);
        void                                setSiteIndex(size_t i);
        void                                print(void) const ;
        
        bool                                operator<(const CharacterEvent& lhs) const;
        bool                                operator>(const CharacterEvent& lhs) const;
        
    protected:
        CharacterEvent(void);
        CharacterEvent(size_t ch_ind, double a, size_t t=CharacterEventType::UNDEFINED);

        size_t                              site_index;                                        // from 0 to N character ("sites")
        double                              age;
        size_t                              event_type;
        bool                                missing;
        
    };
}


#endif

