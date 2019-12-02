#include "PopGenInfinitesSites.h"


using namespace RevBayesCore;



PopGenInfiniteSites::PopGenInfiniteSites(const TypedDagNode<Tree> *tree, const TypedDagNode< RbVector<Tree> > *trees, const TypedDagNode<double> *rate, const TypedDagNode< RbVector<double> > *rates, size_t ns ) :
    TypedDistribution<AbstractHomologousDiscreteCharacterData>( NULL ),
    homogeneous_mutation_rates( rate ),
    heterogeneous_mutation_rates( rates ),
    homogeneous_tau( tree ),
    heterogeneous_tau( trees )
{
    
    redrawValue();
}


PopGenInfiniteSites::~PopGenInfiniteSites(void)
{
    
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PopGenInfiniteSites* PopGenInfiniteSites::clone(void) const
{
    
    return new PopGenInfiniteSites( *this );
}


void PopGenInfiniteSites::compress(void)
{
    
}


double PopGenInfiniteSites::computeLnProbability(void)
{
    
    return 0.0;
}


void PopGenInfiniteSites::redrawValue( void )
{
    
}


void PopGenInfiniteSites::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
}
