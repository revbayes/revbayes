#include "TreeChangeEventHandler.h"
#include "TreeChangeEventListener.h"

#include <iostream>

using namespace RevBayesCore;

TreeChangeEventHandler::TreeChangeEventHandler(void) :
    listeners()
{

}

TreeChangeEventHandler::TreeChangeEventHandler(const TreeChangeEventHandler &h) :
    listeners()
{
    
}


TreeChangeEventHandler::~TreeChangeEventHandler(void) 
{
    
//    if ( listeners.empty() == false )
//    {
//        std::cerr << "Deleting handler while " << listeners.size() << " listener are still listening." << std::endl;
//    }
    
}


TreeChangeEventHandler& TreeChangeEventHandler::operator=(const TreeChangeEventHandler &h)
{
    
    if ( this != &h )
    {
        listeners.clear();
    }
    
    return *this;
}


void TreeChangeEventHandler::addListener(TreeChangeEventListener *l) 
{
    
    listeners.insert( l );
    
}



void TreeChangeEventHandler::fire(const TopologyNode &n, const unsigned& m)
{

    for (std::set<TreeChangeEventListener*>::iterator it = listeners.begin(); it != listeners.end(); ++it) 
    {
        TreeChangeEventListener *l = *it;
        l->fireTreeChangeEvent( n, m );
    }
    
}


const std::set<TreeChangeEventListener*>& TreeChangeEventHandler::getListeners( void ) const 
{
    return listeners;
}


bool TreeChangeEventHandler::isListening(TreeChangeEventListener *l) const
{
    
    // search the set of listeners
    std::set<TreeChangeEventListener*>::iterator pos = listeners.find( l );
    
    return pos != listeners.end();
}


void TreeChangeEventHandler::removeListener(TreeChangeEventListener *l) 
{
    std::set<TreeChangeEventListener*>::iterator pos = listeners.find( l );
    if ( pos != listeners.end() ) 
    {
        listeners.erase( l );
    }
}
