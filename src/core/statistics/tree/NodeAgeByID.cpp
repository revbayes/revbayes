#include "NodeAgeByID.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;

NodeAgeByID::NodeAgeByID(const TypedDagNode<Tree> *t, const size_t qID, const bool s) : TypedFunction<double>( new double(0.0) ),
    tree( t ),
    stemAge( s ),
	queryID( qID )
{

    // add the tree parameter as a parent
    addParameter( tree );
    
    initialize();
    update();
}


NodeAgeByID::~NodeAgeByID( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!

}



NodeAgeByID* NodeAgeByID::clone( void ) const
{
    
    return new NodeAgeByID( *this );
}


void NodeAgeByID::initialize( void )
{
    
	
//	std::cout << queryID << std::endl;
	
}


void NodeAgeByID::update( void )
{
	
    
    if ( stemAge && queryID != tree->getValue().getRoot().getIndex() )
    {
        size_t parentIndex = tree->getValue().getNode(queryID).getParent().getIndex();
        double tmrca = tree->getValue().getNode( parentIndex ).getAge();
        *value = tmrca;
    }
    else
    {
        double tmrca = tree->getValue().getNode( queryID ).getAge();
        *value = tmrca;
    }
}



void NodeAgeByID::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tree) 
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    
}


