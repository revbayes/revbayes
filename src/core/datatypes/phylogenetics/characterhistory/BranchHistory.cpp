#include "BranchHistory.h"

#include <iostream>
#include <iomanip>

#include "CharacterEvent.h"
#include "CharacterEventCompare.h"

using namespace RevBayesCore;


//BranchHistory::BranchHistory(void) : n_characters(0), n_states(0), isTip(false), isRoot(false) { }

BranchHistory::BranchHistory(size_t nc, size_t idx) :
    n_characters(nc),
    branch_index(idx)
{

    parent_characters.resize(n_characters);
    child_characters.resize(n_characters);

}

BranchHistory::BranchHistory(size_t nc, size_t idx, std::set<int> sc) :
    n_characters(nc),
    branch_index(idx)
{

    parent_characters.resize(n_characters);
    child_characters.resize(n_characters);

}



BranchHistory::BranchHistory(const BranchHistory& h) :
    n_characters( h.n_characters ),
    parent_characters(  ),
    child_characters(  ),
    history( h.history ),
    branch_index( h.branch_index)
{

    for (size_t i=0; i<h.parent_characters.size(); ++i)
    {
        CharacterEvent *ce = h.parent_characters[i];
        parent_characters.push_back( ce->clone() );
    }
    
    for (size_t i=0; i<h.child_characters.size(); ++i)
    {
        CharacterEvent *ce = h.child_characters[i];
        child_characters.push_back( ce->clone() );
    }
    
}

BranchHistory::~BranchHistory(void)
{

    for (size_t i = 0; i < parent_characters.size(); i++)
        delete parent_characters[i];
    for (size_t i = 0; i < child_characters.size(); i++)
        delete child_characters[i];
    
//    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it;

}

BranchHistory& BranchHistory::operator=(const BranchHistory &bh)
{

    if (this != &bh)
    {

        for (size_t i = 0; i < parent_characters.size(); i++)
            delete parent_characters[i];
        for (size_t i = 0; i < child_characters.size(); i++)
            delete child_characters[i];
        
        parent_characters.clear();
        child_characters.clear();
        
        
        for (size_t i=0; i<bh.parent_characters.size(); ++i)
        {
            CharacterEvent *ce = bh.parent_characters[i];
            parent_characters.push_back( ce->clone() );
        }
        
        for (size_t i=0; i<bh.child_characters.size(); ++i)
        {
            CharacterEvent *ce = bh.child_characters[i];
            child_characters.push_back( ce->clone() );
        }
        
        n_characters            = bh.n_characters;
        history                 = bh.history;
        branch_index            = bh.branch_index;
    }

    return *this;
}



bool BranchHistory::operator<(const BranchHistory& m) const
{
    return (this < &m);
}


void BranchHistory::addEvent(CharacterEvent* evt)
{
    history.insert(evt);
}

bool BranchHistory::areEventTimesValid(const Tree& tree, const TopologyNode &node) const
{
    // MJL: This causes a problem for branches subtending the root node
    //      if the Newick string does not define a root node branch!
    double lower_boundary = node.getAge();
    double upper_boundary = lower_boundary + tree.getBranchLengthForNode(node);
    double event_age;
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it;
    for (it = history.begin(); it != history.end(); it++)
    {
        event_age = (*it)->getAge();
        if ( event_age > upper_boundary || event_age < lower_boundary)
        {
            return false;
        }
    }
    
    return true;
}


const std::vector<CharacterEvent*>& BranchHistory::getChildCharacters(void) const
{
    return child_characters;
}

std::vector<CharacterEvent*>& BranchHistory::getChildCharacters(void)
{
    return child_characters;
}


std::multiset<CharacterEvent*,CharacterEventCompare>& BranchHistory::getHistory(void)
{
    return history;
}

const std::multiset<CharacterEvent*,CharacterEventCompare>& BranchHistory::getHistory(void) const
{
    return history;
}

const size_t BranchHistory::getNumberCharacters(void) const
{
    return n_characters;
}


const size_t BranchHistory::getNumberEvents(void) const
{
    return history.size();
}


std::vector<CharacterEvent*>& BranchHistory::getParentCharacters(void)
{
    return parent_characters;
}


const std::vector<CharacterEvent*>& BranchHistory::getParentCharacters(void) const
{
    return parent_characters;
}

void BranchHistory::clearEvents(void)
{
    history.clear();
}

void BranchHistory::clearEvents(const std::set<size_t>& indexSet)
{
    
    std::set<size_t>::iterator it_idx;
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h, it_tmp;
    std::set<CharacterEvent*> to_be_deleted;
    
    // for each event in history, delete if index matches indexSet
    for (it_h = history.begin(); it_h != history.end(); it_h++)
    {
        if ( indexSet.find( (*it_h)->getSiteIndex() ) != indexSet.end() )
        {
            to_be_deleted.insert(*it_h);
        }
    }
    
    std::set<CharacterEvent*>::iterator it_d;
    for (it_d = to_be_deleted.begin(); it_d != to_be_deleted.end(); it_d++)
    {
        history.erase(*it_d);
    }

//            it_tmp = it_h;
//            ++it_tmp;
//            history.erase(it_h);
//            //delete *it_h;
//            it_h = it_tmp;
//        }
//        else
//        {
//            ++it_h;
//        }
//    }
}

void BranchHistory::removeEvent(CharacterEvent* evt)
{

    history.erase(evt);

}

void BranchHistory::updateHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet, const std::set<CharacterEvent*>& parentSet, const std::set<CharacterEvent*>& childSet, const std::set<size_t>& indexSet)
{

    /*
    // erase events on branchHistory for indices in indexSet
    clearEvents(indexSet);

    // insert elements into history
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    for (it_h = updateSet.begin(); it_h != updateSet.end(); it_h++)
        history.insert(*it_h);
    */

    updateHistory(updateSet, indexSet);

    // update events on terminal vectors
    std::set<CharacterEvent*>::iterator it_idx;
    for (it_idx = parentSet.begin(); it_idx != parentSet.end(); it_idx++)
        parent_characters[ (*it_idx)->getSiteIndex() ] = *it_idx;
    for (it_idx = childSet.begin(); it_idx != childSet.end(); it_idx++)
        child_characters[ (*it_idx)->getSiteIndex() ] = *it_idx;

}

void BranchHistory::updateHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet, const std::set<size_t>& indexSet)
{
    // erase events on branchHistory for indices in indexSet
    clearEvents(indexSet);

    // insert elements into history
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    for (it_h = updateSet.begin(); it_h != updateSet.end(); it_h++)
    {
        history.insert(*it_h);
    }

}

void BranchHistory::updateHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& updateSet)
{
    // erase events on branchHistory for indices in indexSet
    clearEvents();

    // insert elements into history
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    for (it_h = updateSet.begin(); it_h != updateSet.end(); it_h++)
    {
        history.insert(*it_h);
    }

}

void BranchHistory::setChildCharacters(const std::vector<CharacterEvent*>& s)
{
    std::vector<CharacterEvent*>::const_iterator it;
    for (it = s.begin(); it != s.end(); it++)
    {
        child_characters[ (*it)->getSiteIndex() ] = *it;//new CharacterEvent(**it);
    }

}

void BranchHistory::setChildCharacters(const std::set<CharacterEvent*>& s)
{
    std::set<CharacterEvent*>::iterator it;
    for (it = s.begin(); it != s.end(); it++)
    {
        child_characters[ (*it)->getSiteIndex() ] = *it;//new CharacterEvent(**it);
    }

}

void BranchHistory::setParentCharacters(const std::vector<CharacterEvent*>& s)
{
    std::vector<CharacterEvent*>::const_iterator it;
    for (it = s.begin(); it != s.end(); it++)
    {
        parent_characters[ (*it)->getSiteIndex() ] = *it;//new CharacterEvent(**it);
    }

}

void BranchHistory::setParentCharacters(const std::set<CharacterEvent*>& s)
{
    std::set<CharacterEvent*>::iterator it;
    for (it = s.begin(); it != s.end(); it++)
    {
        parent_characters[ (*it)->getSiteIndex() ] = *it;//new CharacterEvent(**it);
    }

}


void BranchHistory::setHistory(const std::set<CharacterEvent*,CharacterEventCompare>& s)
{
    history.clear();
    for (std::set<CharacterEvent*,CharacterEventCompare>::iterator it = s.begin(); it != s.end(); it++)
    {
        history.insert(*it);
    }

}

void BranchHistory::setHistory(const std::multiset<CharacterEvent*,CharacterEventCompare>& s)
{
    history = s;
}


void BranchHistory::print(const TopologyNode* nd) const
{
    
    std::set<CharacterEvent*,CharacterEventCompare>::iterator it_h;
    std::vector<CharacterEvent*>::iterator it_v;

    double start_age=0.0;
    double end_age=0.0;
    if (nd != NULL) {
        end_age = nd->getAge();
        if (!nd->isRoot()) {
            start_age = nd->getParent().getAge();
        } else {
            start_age = 1e6;
        }
    }
//    std::cout << parentCharacters.size() << "\n";
//    for (size_t i = 0; i < parentCharacters.size(); i++)
//        std::cout << parentCharacters[i] << " ";
//    std::cout << "\n";

    std::vector<CharacterEvent*> tmp = parent_characters;

    std::cout << "BranchHistory " << branch_index << " size=" << history.size() << "  " << this << "\n";
    std::cout << "                             ";
    for (size_t i = 0; i < n_characters; i++)
    {
        if (i % 10 == 0) std::cout << ".";
        else std::cout << " ";
    }
    std::cout << "\n";
    std::cout << "end           ";
    std::cout << std::setw(12) << std::setprecision(6) << end_age << " : ";
//    std::cout << "                       end : ";
    for (it_v = child_characters.begin(); it_v != child_characters.end(); it_v++)
    {
//        std::cout << (*it_v)->getState();
        std::cout << (*it_v)->getStateStr();
    }
    std::cout << "\n";

    for (it_h = history.begin(); it_h != history.end(); it_h++)
    {
        std::cout << *it_h << "   ";
        std::cout << std::setw(12) << std::setprecision(6) << (*it_h)->getAge() << " : ";
        tmp[ (*it_h)->getSiteIndex() ] = *it_h;
        for (size_t i = 0; i < n_characters; i++)
        {
            if (i != (*it_h)->getSiteIndex())
            {
                std::cout << " ";
            }
            else
            {
                std::cout << (*it_h)->getStateStr();
            }

//                std::cout << (*it_h)->getState();
            //std::cout << " ";
        }
        std::cout << "\n";

    }
//    std::cout << "                     start : ";
    std::cout << "start         ";
    std::cout << std::setw(12) << std::setprecision(6) << start_age << " : ";
    for (it_v = parent_characters.begin(); it_v != parent_characters.end(); it_v++)
    {
//        std::cout << (*it_v)->getState();
        std::cout << (*it_v)->getStateStr();
    }
    std::cout << "\n";
    std::cout << "                             ";
    for (size_t i = 0; i < n_characters; i++)
    {
        if (i % 10 == 0)
            std::cout << ".";
        else
            std::cout << " ";
    }
    std::cout << "\n";
    ;
}


const size_t BranchHistory::getBranchIndex(void) const
{
    return branch_index;
}


CharacterEvent* BranchHistory::getEvent(size_t i)
{
    std::multiset<CharacterEvent*,CharacterEventCompare>::iterator it = history.begin();
    for (size_t j=0; j<i; ++j)
    {
        ++it;
    }

    return (*it);
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const BranchHistory& x)
{

    o << x.getNumberEvents();

    return o;
}
