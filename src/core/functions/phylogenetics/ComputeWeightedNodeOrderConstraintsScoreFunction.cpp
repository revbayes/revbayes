#include "ComputeWeightedNodeOrderConstraintsScoreFunction.h"

#include <cstddef>
#include <iostream>
#include <vector>

#include "RbException.h"
#include "TreeUtilities.h"
#include "Cloneable.h"
#include "RelativeNodeAgeWeightedConstraints.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

computeWeightedNodeOrderConstraintsScoreFunction::computeWeightedNodeOrderConstraintsScoreFunction(const TypedDagNode<Tree> *t, const TypedDagNode<RelativeNodeAgeWeightedConstraints> *co, TypedDagNode< double>* bet) : TypedFunction<double>( new double(0.0) ),
tau( t ),
constraints( co ),
beta( bet ),
nodeAges(),
constrainedNodes()
{
    try {
    // add the constraints parameter as a parent
    addParameter( tau );
    addParameter( constraints );
    addParameter( beta );
    updateSetOfConstrainedNodes();

    update();
}
catch (RbException &e)
{
    std::cerr << e.getMessage() << std::endl;
    }

}


computeWeightedNodeOrderConstraintsScoreFunction::~computeWeightedNodeOrderConstraintsScoreFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



computeWeightedNodeOrderConstraintsScoreFunction* computeWeightedNodeOrderConstraintsScoreFunction::clone( void ) const
{
    return new computeWeightedNodeOrderConstraintsScoreFunction( *this );
}


void computeWeightedNodeOrderConstraintsScoreFunction::keep( const DagNode *affecter)
{
    //delegate to base class
    TypedFunction< double >::keep( affecter );

}


void computeWeightedNodeOrderConstraintsScoreFunction::reInitialized( void )
{
    *value = 0.0;
}


void computeWeightedNodeOrderConstraintsScoreFunction::restore( const DagNode *restorer)
{
    //delegate to base class
    TypedFunction< double >::restore( restorer );
}


void computeWeightedNodeOrderConstraintsScoreFunction::touch(const DagNode *toucher)
{

    //delegate to base class
    TypedFunction< double >::touch( toucher );

}


double computeWeightedNodeOrderConstraintsScoreFunction::computeWeightedNodeAgeConstraintsScore ( void )
{
    std::vector < std::pair < std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double > > constra = constraints->getValue().getConstraints();

    double score = 0.0;
    for (size_t i = 0; i < constra.size() ; ++i) {
        if ( nodeAges.at(constra[i].first.first) <  nodeAges.at(constra[i].first.second) ) {
            score += beta->getValue() * constra[i].second * ( nodeAges.at(constra[i].first.first) - nodeAges.at(constra[i].first.second) ) ;
        }
    }
    return score;

}


void computeWeightedNodeOrderConstraintsScoreFunction::update( void )
{
    (*value) = 0.0;
    try {
        updateMapOfNodeAges();
    }
    catch (RbException &e)
    {
        std::cerr << e.getMessage() << std::endl;
    }

    (*value) = computeWeightedNodeAgeConstraintsScore() ;

    return;

}


void computeWeightedNodeOrderConstraintsScoreFunction::updateSetOfConstrainedNodes()
{
    std::vector < std::pair < std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> >, double > > constra = constraints->getValue().getConstraints();
    for (size_t i = 0; i < constra.size() ; ++i) {
        constrainedNodes.insert(constra[i].first.first);
        constrainedNodes.insert(constra[i].first.second);
    }
    return;
}




//Here we compute node ages from the current tree.
void computeWeightedNodeOrderConstraintsScoreFunction::updateMapOfNodeAges()
{


    nodeAges.clear();
    for (std::set< std::pair < std::string, std::string > >::iterator elem=constrainedNodes.begin(); elem != constrainedNodes.end(); ++elem)
    {
        nodeAges[(*elem)] = TreeUtilities::getAgeOfMRCA(tau->getValue(), elem->first, elem->second);
    }

    return;

}


void computeWeightedNodeOrderConstraintsScoreFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == constraints)
    {
        constraints = static_cast<const TypedDagNode<RelativeNodeAgeWeightedConstraints>* >( newP );
        updateSetOfConstrainedNodes();
    }
}
