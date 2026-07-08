#include <cstddef>
#include <algorithm>
#include <iosfwd>
#include <string>
#include <vector>

#include "DagNode.h"
#include "Distribution.h"
#include "RbException.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "Parallelizable.h"

namespace RevBayesCore { template <class valueType> class RbOrderedSet; }

using namespace RevBayesCore;



Distribution::Distribution(void) : Parallelizable(),
    parameters()
{
    
}


Distribution::Distribution(const Distribution &d) : Parallelizable(d),
    parameters( d.parameters )
{
    
    for (std::vector<const DagNode*>::iterator it=parameters.begin(); it!=parameters.end(); ++it)
    {
        (*it)->incrementReferenceCount();
    }
    
}


Distribution::~Distribution( void )
{
    
    for (std::vector<const DagNode*>::iterator it=parameters.begin(); it!=parameters.end(); ++it)
    {
        const DagNode *the_node = *it;
        if ( the_node->decrementReferenceCount() == 0 )
        {
            delete the_node;
        }
    }
    
}



Distribution& Distribution::operator=(const Distribution &d)
{
    Parallelizable::operator=(d);
    
    if ( this != &d )
    {
        
        for (std::vector<const DagNode*>::iterator it=parameters.begin(); it!=parameters.end(); ++it)
        {
            const DagNode *the_node = *it;
            if ( the_node->decrementReferenceCount() == 0 )
            {
                delete the_node;
            }
        }
        parameters.clear();
        
        parameters = d.parameters;
        
        for (std::vector<const DagNode*>::iterator it=parameters.begin(); it!=parameters.end(); ++it)
        {
            (*it)->incrementReferenceCount();
        }
    }
    
    return *this;
}


/**
 * Method stub: override for specialized treatment. 
 */
void Distribution::bootstrap( void )
{
    // do nothing
}

/*
 * By default, distribution parameter changes do not affect child values.
 */
bool Distribution::childrenAreAffectedBy(const DagNode *) const
{
    return false;
}


/**
 * Add this parameter to our set of parameters.
 */
void Distribution::addParameter(const DagNode *p)
{
    
    // only if the parameter is not NULL
    if ( p != NULL )
    {
        std::vector<const DagNode*>::iterator pos = std::find(parameters.begin(), parameters.end(), p);
        if ( pos == parameters.end() )
        {
            parameters.push_back( p );
            
            // increment reference count
            p->incrementReferenceCount();
            
        }
    
    }
    
}



RevLanguage::RevPtr<RevLanguage::RevVariable> Distribution::executeProcedure(const std::string &n, const std::vector<DagNode *> args, bool &f)
{
    // no function found
    return NULL;
}


std::vector<double> Distribution::getMixtureProbabilities(void) const
{
    return std::vector<double>(1, 1.0);
}


size_t Distribution::getNumberOfMixtureElements(void) const
{
    return 0;
}


/**
 * Get a const reference to the set of parameters for this distribution.
 */
const std::vector<const DagNode*>& Distribution::getParameters( void ) const
{
    
    return parameters;
}


/* Invalidate local cached state; fullyInvalidateSelf disables partial invalidation. */
void Distribution::invalidate( const DagNode *toucher, bool fullyInvalidateSelf )
{
    invalidateSpecialization( toucher, fullyInvalidateSelf );
}


/* Method stub: override for specialized local cache invalidation. */
void Distribution::invalidateSpecialization( const DagNode *, bool )
{
    // do nothing
}


/* Method stub: override for specialized treatment. */
void Distribution::keep( void )
{
    // delegate to specialization
    keepSpecialization();
    
}


/* Method stub: override for specialized treatment. */
void Distribution::keepSpecialization(void)
{
    // do nothing
}


/* Method stub: override for specialized treatment. */
void Distribution::reInitialized( void )
{
    // do nothing
}


/**
 * Remove this parameter from our list of parameters.
 */
void Distribution::removeParameter(const RevBayesCore::DagNode *p)
{
    
    // only if the parameter is not NULL
    if ( p != NULL )
    {
     
        std::vector<const DagNode *>::iterator it = std::find( parameters.begin(), parameters.end(), p );
        if ( it != parameters.end() )
        {
            parameters.erase( it );
            if ( p->decrementReferenceCount() == 0 )
            {
                delete p;
            }
            
        }
    
    }
    
}


/* Method stub: override for specialized treatment. */
void Distribution::restore( void )
{
    // delegate to specialization
    restoreSpecialization();
}


/* Method stub: override for specialized treatment. */
void Distribution::restoreSpecialization(void)
{
    // do nothing
}


/**
 * Set the flag if we are in MCMC mode or not.
 * This default implementation does nothing but derived classes can simply overwrite this.
 * It is just a placeholder.
 */
void Distribution::setMcmcMode(bool tf)
{
    // do nothing
}


/* Method stub: override for specialized treatment. */
void Distribution::snapshot( void )
{
    snapshotSpecialization();
}


/* Method stub: override for specialized treatment. */
void Distribution::snapshotSpecialization( void )
{
    // do nothing
}


/**
 * Swap the old parameter with a new one.
 * This will be called for example when the entire model graph is cloned or
 * when we replace a variable with the same name (re-assignment).
 * Here we update our set and delegate to the derived class.
 */
void Distribution::swapParameter(const DagNode *old_p, const DagNode *new_p)
{
    
    std::vector<const DagNode *>::iterator position = std::find(parameters.begin(), parameters.end(), old_p);
    if ( position != parameters.end() )
    {
//        parameters.erase( position );
//        parameters.push_back( newP );
        (*position) = new_p;
        swapParameterInternal( old_p, new_p );
        
        // increment and decrement the reference counts
        new_p->incrementReferenceCount();
        if ( old_p->decrementReferenceCount() == 0 )
        {
            throw RbException("Memory leak in distribution. Please report this bug to Sebastian.");
        }
        
    }
    else
    {
        
        throw RbException() << "Could not find the distribution parameter to be swapped: " << old_p->getName() << " to " << new_p->getName(); 
    
    }
    
}
