#include "TreeAssemblyFunction.h"

#include <cstddef>
#include <vector>

#include "RbException.h"
#include "DagNode.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

TreeAssemblyFunction::TreeAssemblyFunction(const TypedDagNode<Tree> *t, const TypedDagNode< RbVector<double> > *b) : TypedFunction<Tree>( NULL ),
    tau( t ),
    brlen( b )
{

    if (tau->getValue().getNumberOfNodes() - 1 != brlen->getValue().size())
    {
        throw(RbException("Number of branches does not match the number of branch lengths"));
    }

    // add the lambda parameter as a parent
    addParameter( tau );
    addParameter( brlen );
    
    value = const_cast<Tree*>( &tau->getValue() );
    
    brlenFlagDirty = true;
    update();
}


TreeAssemblyFunction::TreeAssemblyFunction(const TreeAssemblyFunction &f) : TypedFunction<Tree>( f ),
    tau( f.tau ),
    brlen( f.brlen )
{
    
    // the base class has created a new value instance
    // we need to delete it here to avoid memory leaks
    delete value;
    
    value = const_cast<Tree*>( &tau->getValue() );
    
    brlenFlagDirty = true;
    update();
}


TreeAssemblyFunction::~TreeAssemblyFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
    
    // rescue deletion
    value = NULL;
}



TreeAssemblyFunction* TreeAssemblyFunction::clone( void ) const
{
    return new TreeAssemblyFunction( *this );
}


void TreeAssemblyFunction::keep( const DagNode *affecter )
{
    //delegate to base class
    TypedFunction< Tree >::keep( affecter );
    
    // SH (20200221): This needs to stay, otherwise the MCMC does not update properly.
    touchedNodeIndices.clear();
    brlenFlagDirty = false;
    
    // SH (20200221): We must not call update because it breaks the MCMC!
    // SH (20190913): There seems to be an issue if we use two replicates
    // So we need to make sure keep is only called after update!
//    update();
    
}


void TreeAssemblyFunction::reInitialized( void )
{
    
}


void TreeAssemblyFunction::restore( const DagNode *restorer )
{
    //delegate to base class
    TypedFunction< Tree >::restore( restorer );
    
    touchedNodeIndices.clear();
    brlenFlagDirty = false;
}


void TreeAssemblyFunction::touch(const DagNode *toucher)
{
    
    //delegate to base class
    TypedFunction< Tree >::touch( toucher );
    
    //reset flag
    brlenFlagDirty = true;
    
    if ( toucher == brlen )
    {
        const std::set<size_t> &touchedIndices = toucher->getTouchedElementIndices();
        touchedNodeIndices.insert(touchedIndices.begin(), touchedIndices.end());
    }
    
    if ( toucher != tau && touchedNodeIndices.size() == 0 )
    {
        for (size_t i = 0; i < brlen->getValue().size(); ++i)
        {
            touchedNodeIndices.insert(i);
        }
    }
}


void TreeAssemblyFunction::update( void )
{

    const std::vector<double> &v = brlen->getValue();
    if ( touchedNodeIndices.size() < v.size() )
    {
        for (size_t i = 0; i < v.size(); ++i)
        {
            value->getNode(i).setBranchLength( v[i], false );
        }
    }

    if ( touchedNodeIndices.size() > 0 )
    {
        for (std::set<size_t>::iterator it = touchedNodeIndices.begin(); it != touchedNodeIndices.end(); ++it)
        {
            value->getNode(*it).setBranchLength( v[*it], brlenFlagDirty );
        }
        touchedNodeIndices.clear();
        brlenFlagDirty = false;
    }

}



void TreeAssemblyFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
        
        Tree *psi = const_cast<Tree*>( &tau->getValue() );
        
        // finally store the new value
        value = psi;
        
    }
    else if (oldP == brlen)
    {
        brlen = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}


