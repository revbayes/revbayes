#include "TreeAssemblyFunction.h"
#include "RbException.h"

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


void TreeAssemblyFunction::keep(DagNode *affecter)
{
    //delegate to base class
    TypedFunction< Tree >::keep( affecter );
    
//    touchedNodeIndices.clear();
}


void TreeAssemblyFunction::reInitialized( void )
{
    
}


void TreeAssemblyFunction::restore(DagNode *restorer)
{
    //delegate to base class
    TypedFunction< Tree >::restore( restorer );
    
//    touchedNodeIndices.clear();
}


void TreeAssemblyFunction::touch(DagNode *toucher)
{
    
    //delegate to base class
    TypedFunction< Tree >::touch( toucher );
    
    if ( toucher == brlen )
    {
        const std::set<size_t> &touchedIndices = toucher->getTouchedElementIndices();
        touchedNodeIndices.insert(touchedIndices.begin(), touchedIndices.end());
    }
    
}


void TreeAssemblyFunction::update( void )
{
    
    if ( touchedNodeIndices.size() > 0 )
    {
        const std::vector<double> &v = brlen->getValue();
        for (std::set<size_t>::iterator it = touchedNodeIndices.begin(); it != touchedNodeIndices.end(); ++it)
        {
            value->getNode(*it).setBranchLength(v[*it]);
        }
        touchedNodeIndices.clear();
    }
    else
    {
        const std::vector<double> &v = brlen->getValue();
        for (size_t i = 0; i < v.size(); ++i)
        {
            value->getNode(i).setBranchLength( v[i] );
        }
        
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


