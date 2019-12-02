#include "BiogeographyRateGeneratorSequenceFunction.h"
#include "ConstantNode.h"
#include "RateGenerator.h"
#include "RateMatrix_JC.h"
#include "RbException.h"

using namespace RevBayesCore;

BiogeographyRateGeneratorSequenceFunction::BiogeographyRateGeneratorSequenceFunction(size_t nc, bool fe, unsigned mrs) : TypedFunction<RateGeneratorSequence>( new RateGeneratorSequence_Biogeography( nc, fe, mrs ) )
{
//    homogeneousGainLossRates            = new ConstantNode<RbVector<double> >("homogeneousGainLossRates", new RbVector<double>(2,0.5));
//    heterogeneousGainLossRates          = NULL;
    homogeneousRateMatrix               = new ConstantNode<RateGenerator>("homogeneousRateMatrix", new RateMatrix_JC(2));
    
    useGeographicDistance               = false;
    
    this->addParameter(homogeneousRateMatrix);
    
    update();
}


BiogeographyRateGeneratorSequenceFunction::~BiogeographyRateGeneratorSequenceFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



BiogeographyRateGeneratorSequenceFunction* BiogeographyRateGeneratorSequenceFunction::clone( void ) const {
    return new BiogeographyRateGeneratorSequenceFunction( *this );
}


void BiogeographyRateGeneratorSequenceFunction::update( void )
{
    
    // touch specialization for granular updates?
    
    // set the gainLossRate
    const RateGenerator& rm = homogeneousRateMatrix->getValue();
    static_cast< RateGeneratorSequence_Biogeography* >(value)->setRateMatrix(&rm);

    
    // set the distancePower
    if (useGeographicDistance)
    {
        ;
//        const DistanceRateModifier& drm = const_cast<const DistanceRateModifier>( distanceRateModifier->getValue() );
//        static_cast< RateGeneratorSequence_Biogeography* >(value)->setDistanceRateModifier(drm);
    }

    // @Michael: Might add back root frequencies!
//    if (rootFrequencies != NULL)
//    {
//        const std::vector<double>& f = rootFrequencies->getValue();
//        static_cast<RateGeneratorSequence_Biogeography*>(value)->setRootFrequencies(f);
//    }
    
    value->updateMap();
}

void BiogeographyRateGeneratorSequenceFunction::setRateMatrix(const TypedDagNode<RateGenerator>* r)
{
    // remove the old parameter first
    if ( homogeneousRateMatrix != NULL )
    {
        this->removeParameter( homogeneousRateMatrix );
        homogeneousRateMatrix = NULL;
    }
    
    // set the value
    homogeneousRateMatrix = r;
    
    // add the parameter
    this->addParameter( homogeneousRateMatrix );
}


void BiogeographyRateGeneratorSequenceFunction::setDistanceRateModifier(const TypedDagNode<CharacterHistoryRateModifier> *drm)
{
    
    // remove the old parameter first
    if ( distanceRateModifier != NULL )
    {
        this->removeParameter( distanceRateModifier );
        distanceRateModifier = NULL;
    }
    
    // set the value
    useGeographicDistance = true;
    distanceRateModifier = drm;
    
    // add the parameter
    this->addParameter( distanceRateModifier );

}


void BiogeographyRateGeneratorSequenceFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == homogeneousRateMatrix)
    {
        homogeneousRateMatrix = static_cast<const TypedDagNode<RateGenerator>* >( newP );
    }
    else if (oldP == distanceRateModifier)
    {
        distanceRateModifier = static_cast<const TypedDagNode<CharacterHistoryRateModifier>* >( newP );
    }

}
