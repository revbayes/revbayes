#include "CharacterEvent.h"
#include "CharacterEventType.h"
#include <iostream>
#include <sstream>

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

CharacterEvent::CharacterEvent(const CharacterEvent& c) :
    site_index(c.site_index),
    age(c.age),
    event_type(c.event_type)
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


void CharacterEvent::setEventType(size_t t)
{
    event_type = t;
}


void CharacterEvent::setSiteIndex(size_t i)
{
    site_index = i;
}


void CharacterEvent::print(void) const
{
    std::cout << site_index << " " << getStateStr() << " " << age << "\n";
}

