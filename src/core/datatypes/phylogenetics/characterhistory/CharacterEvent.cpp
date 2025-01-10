#include "CharacterEvent.h"

#include <iostream>

using namespace RevBayesCore;


CharacterEvent::CharacterEvent(void)
{
    
}

CharacterEvent::CharacterEvent(size_t ch_ind, double a, size_t t) :
    site_index(ch_ind),
    age(a),
    event_type(t)
{
    
}

CharacterEvent::~CharacterEvent(void)
{
    
}

bool CharacterEvent::operator<(const CharacterEvent& rhs) const
{
    return age < rhs.age;
}

bool CharacterEvent::operator>(const CharacterEvent& rhs) const
{
    return age > rhs.age;
}

double CharacterEvent::getAge(void) const
{
    return age;
}

size_t CharacterEvent::getEventType(void) const
{
    return event_type;
}

size_t CharacterEvent::getSiteIndex(void) const
{
    return site_index;
}


void CharacterEvent::setAge(double a)
{
    age = a;
}


void CharacterEvent::setSiteIndex(size_t i)
{
    site_index = i;
}


void CharacterEvent::setMissingState(bool tf)
{
    missing = tf;
}


void CharacterEvent::setEventType(size_t t)
{
    event_type = t;
}


void CharacterEvent::print(void) const
{
    std::cout << site_index << " " << getStateStr() << " " << age << "\n";
}

