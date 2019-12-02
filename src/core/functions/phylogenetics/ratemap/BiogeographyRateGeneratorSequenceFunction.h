#ifndef BiogeographyRateGeneratorSequenceFunction_H
#define BiogeographyRateGeneratorSequenceFunction_H

//#include "CountRateModifier.h"
#include "DistanceRateModifier.h"
#include "RateGenerator.h"
#include "RateGeneratorSequence_Biogeography.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    class BiogeographyRateGeneratorSequenceFunction : public TypedFunction<RateGeneratorSequence> {
        
    public:
        BiogeographyRateGeneratorSequenceFunction(size_t nc, bool fe=true, unsigned mrs=0);                                                                                             // pass in geography object??
        virtual                                            ~BiogeographyRateGeneratorSequenceFunction(void);                                              //!< Virtual destructor
        
        // public member functions
        BiogeographyRateGeneratorSequenceFunction*                        clone(void) const;                                                              //!< Create an independent clone
        
        // set parameters
        void                                                setDistanceRateModifier(const TypedDagNode<CharacterHistoryRateModifier>* drm);
//        void                                                setCountRateModifier(const TypedDagNode<double>* rsrm);
        void                                                setGainLossRates(const TypedDagNode< RbVector<double> >* glr);
        void                                                setRateMatrix(const TypedDagNode<RateGenerator>* rm);
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<RateGenerator>*                  homogeneousRateMatrix;

        const TypedDagNode<CharacterHistoryRateModifier>*   distanceRateModifier;
//      const TypedDagNode<CountRateModifier>*              countRateModifier;
        
        // geography epochs
        // rate epochs
        // branchwise models, etc
        
        bool                                                useGeographicDistance;
        
    };
    
}

#endif
