#include "RootedTripletDistributionFunction.h"
#include "RbException.h"
#include "RootedTripletDistribution.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TreeUtilities.h"

#include <boost/dynamic_bitset.hpp>

using namespace RevBayesCore;

RootedTripletDistributionFunction::RootedTripletDistributionFunction( const TypedDagNode< RbVector< Tree > > *ts,  const TypedDagNode< RbVector< std::string > > *sn ) : TypedFunction<RootedTripletDistribution>( new RootedTripletDistribution(  ) ) {
//    ts->getValue(), sn->getValue()
    trees = ts ;
    species = sn->clone();
    rtd = new RootedTripletDistribution( trees->getValue(), species->getValue()  );
    // add the lambda parameter as a parent
    addParameter( ts );
   // value->setSpecies ( species->getValue() );
    
   // update();
}


RootedTripletDistributionFunction::RootedTripletDistributionFunction( const TypedDagNode< RbVector< std::string > > *sn ): TypedFunction<RootedTripletDistribution>( new RootedTripletDistribution(  ) ) {
    //    ts->getValue(), sn->getValue()
    species = sn->clone();
    rtd = new RootedTripletDistribution( species->getValue()  );

}

RootedTripletDistributionFunction::RootedTripletDistributionFunction( const TypedDagNode<RevBayesCore::RbVector< Taxon > > *t ): TypedFunction<RootedTripletDistribution>( new RootedTripletDistribution(  ) ) {
    //    ts->getValue(), sn->getValue()
    taxa = t->clone();
    rtd = new RootedTripletDistribution( taxa->getValue()  );
    
}


RootedTripletDistributionFunction::~RootedTripletDistributionFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



RootedTripletDistributionFunction* RootedTripletDistributionFunction::clone( void ) const
{
    return new RootedTripletDistributionFunction( *this );
}


void RootedTripletDistributionFunction::update( void )
{
    if (rtd != NULL  ) {
        rtd->resetDistribution();
    }
    
    rtd->setTrees( trees->getValue() );
    
    // TODO: Bastien, this does not compile (Michael/Sebastian)
    //rtd->setRecordBranchLengths(recordBranchLengths->getValue());
    rtd->extractTriplets(  );
    value = rtd->clone();
}



void RootedTripletDistributionFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    trees = static_cast<const TypedDagNode< RbVector<Tree> >* >( newP );

}


void RootedTripletDistributionFunction::setTrees(TypedDagNode< RbVector< Tree > >* ts)
{

    trees = ts;
    removeParameter( ts );
    addParameter( ts );

}


void RootedTripletDistributionFunction::setRecordBranchLengths(  TypedDagNode< bool >* b)
{
    recordBranchLengths = b;
   // update();
}
