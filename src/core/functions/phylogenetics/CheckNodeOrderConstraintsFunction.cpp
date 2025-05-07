#include "CheckNodeOrderConstraintsFunction.h"

#include <cstddef>
#include <iostream>
#include <vector>

#include "RbException.h"
#include "TreeUtilities.h"
#include "Cloneable.h"
#include "RelativeNodeAgeConstraints.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

CheckNodeOrderConstraintsFunction::CheckNodeOrderConstraintsFunction(const TypedDagNode<Tree> *t, const TypedDagNode<RelativeNodeAgeConstraints> *co) : TypedFunction<Boolean>( new Boolean(false) ),
tau( t ),
constraints( co ),
nodeAges(),
constrainedNodes()
{
    try {
    // add the constraints parameter as a parent
    addParameter( tau );
    addParameter( constraints );
    updateSetOfConstrainedNodes();
    
    update();
}
catch (RbException &e)
{
    std::cerr << e.getMessage() << std::endl;
    }

}


CheckNodeOrderConstraintsFunction::~CheckNodeOrderConstraintsFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



CheckNodeOrderConstraintsFunction* CheckNodeOrderConstraintsFunction::clone( void ) const
{
    return new CheckNodeOrderConstraintsFunction( *this );
}


void CheckNodeOrderConstraintsFunction::keep( const DagNode *affecter)
{
    //delegate to base class
    TypedFunction< Boolean >::keep( affecter );
    
}


void CheckNodeOrderConstraintsFunction::reInitialized( void )
{
    *value = Boolean(false);
}


void CheckNodeOrderConstraintsFunction::restore( const DagNode *restorer)
{
    //delegate to base class
    TypedFunction< Boolean >::restore( restorer );
}


void CheckNodeOrderConstraintsFunction::touch(const DagNode *toucher)
{
    
    //delegate to base class
    TypedFunction< Boolean >::touch( toucher );
    
}


bool CheckNodeOrderConstraintsFunction::checkNodeAgeConstraints ( void )
{
    std::vector <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > > constra = constraints->getValue().getConstraints();
    for (size_t i = 0; i < constra.size() ; ++i) {
        if ( nodeAges.at(constra[i].first) <  nodeAges.at(constra[i].second) ) {
            return false;
        }
    }
    return true;
}


void CheckNodeOrderConstraintsFunction::update( void )
{
    (*value) = Boolean(false);
    try {
        updateMapOfNodeAges();
    }
    catch (RbException &e)
    {
        std::cerr << e.getMessage() << std::endl;
    }

    (*value) = Boolean( checkNodeAgeConstraints() );

    return;
    
}


void CheckNodeOrderConstraintsFunction::updateSetOfConstrainedNodes()
{
    std::vector <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > > constra = constraints->getValue().getConstraints();
    for (size_t i = 0; i < constra.size() ; ++i) {
        constrainedNodes.insert(constra[i].first);
        constrainedNodes.insert(constra[i].second);
    }
    return;
}




//Here we compute node ages from the current tree.
void CheckNodeOrderConstraintsFunction::updateMapOfNodeAges()
{
  
    
    nodeAges.clear();
    for (std::set< std::pair < std::string, std::string > >::iterator elem=constrainedNodes.begin(); elem != constrainedNodes.end(); ++elem)
    {
        nodeAges[(*elem)] = TreeUtilities::getAgeOfMRCA(tau->getValue(), elem->first, elem->second);
    }
    
}


void CheckNodeOrderConstraintsFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == constraints)
    {
        constraints = static_cast<const TypedDagNode<RelativeNodeAgeConstraints>* >( newP );
        updateSetOfConstrainedNodes();
    }
}

















