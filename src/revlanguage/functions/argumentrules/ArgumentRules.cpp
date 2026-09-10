/**
 * @file
 * This file contains the declaration of ArgumentRules, which is
 * is the container for all argument rules. We need a special object
 * for proper memory allocation and copying of the objects.
 *
 * @brief Declaration of ArgumentRules
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-02-13, version 1.0
 *
 * $Id$
 */

#include <sstream>
#include "ArgumentRules.h"
#include "ArgumentRule.h"

using namespace RevLanguage;

ArgumentRules::ArgumentRules(void) {
    
}

ArgumentRules::ArgumentRules(const ArgumentRules& ar)
{
    
    // get a copy of the argument rules
    for (std::vector<ArgumentRule*>::const_iterator it = ar.rules.begin(); it != ar.rules.end(); it++)
    {
        rules.push_back( (*it)->clone() );
    }
    
}


ArgumentRules::~ArgumentRules(void)
{
    clear();
}


std::vector<ArgumentRule*>::iterator ArgumentRules::begin(void)
{
    return rules.begin();
}


std::vector<ArgumentRule*>::const_iterator ArgumentRules::begin(void) const
{
    return rules.begin();
}


ArgumentRules& ArgumentRules::operator=(const ArgumentRules &a)
{
    
    if (this != &a)
    {
        // first we clear ourselves
        clear();
        
        // get a copy of the argument rules
        for (std::vector<ArgumentRule*>::const_iterator it = a.rules.begin(); it != a.rules.end(); it++)
        {
            rules.push_back( (*it)->clone() );
        }
    }
    
    return (*this);
}


const ArgumentRule& ArgumentRules::operator[](size_t index) const
{
    return *rules[index];
}


ArgumentRule& ArgumentRules::operator[](size_t index)
{
    return *rules[index];
}


void ArgumentRules::clear(void)
{
    // delete the argument rules
    for (size_t i = 0; i != rules.size(); ++i)
    {
        ArgumentRule* the_rule = rules[i];
        delete the_rule;
    }
    rules.clear();
}

std::vector<ArgumentRule*>::iterator ArgumentRules::end(void)
{
    return rules.end();
}


std::vector<ArgumentRule*>::const_iterator ArgumentRules::end(void) const
{
    return rules.end();
}



void ArgumentRules::insert(std::vector<ArgumentRule *>::iterator it, ArgumentRule *rule)
{
    rules.insert(it, rule);
}


void ArgumentRules::insert(std::vector<ArgumentRule *>::iterator it, std::vector<ArgumentRule *>::const_iterator first, std::vector<ArgumentRule *>::const_iterator last) {
    
    std::vector<ArgumentRule*> tmp;
    for (std::vector<ArgumentRule *>::const_iterator org = first; org != last; org++)
    {
        ArgumentRule* theOriginalRule = *org;
        ArgumentRule* theClone = theOriginalRule->clone();
        tmp.push_back(theClone);
    }
    
    rules.insert(it, tmp.begin(), tmp.end());
}


/** Push back a new rule. We just delegate it to our member variable. */
void ArgumentRules::push_back(ArgumentRule *rule)
{
    rules.push_back(rule);
}


/** Push back a new rule. We just delegate it to our member variable. */
void ArgumentRules::push_front(ArgumentRule *rule)
{
    std::vector<ArgumentRule*>::iterator it = rules.begin();
    rules.insert(it,rule);
}



/** Get the size of the vector. We delegate to the member. */
size_t ArgumentRules::size(void) const
{
    return rules.size();
}


void ArgumentRules::printValue(std::ostream& o) const
{
    o << " (";
    for (size_t i=0; i<size(); i++)
    {
        if (i != 0)
        {
            o << ", ";
        }
        (*this)[i].printValue(o);
    }
    o << ")";
}

std::string ArgumentRules::show() const
{
    std::stringstream o;
    printValue(o);
    return o.str();
}
