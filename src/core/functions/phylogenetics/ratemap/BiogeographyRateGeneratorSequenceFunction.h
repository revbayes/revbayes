#ifndef BiogeographyRateGeneratorSequenceFunction_H
#define BiogeographyRateGeneratorSequenceFunction_H

#include <cstddef>

#include "TypedFunction.h"
#include "RateGeneratorSequence.h"

namespace RevBayesCore {
class CharacterHistoryRateModifier;
class DagNode;
class RateGenerator;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
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
