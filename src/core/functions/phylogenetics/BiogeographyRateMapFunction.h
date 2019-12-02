//#ifndef BiogeographyRateMapFunction_H
//#define BiogeographyRateMapFunction_H
//
////#include "CountRateModifier.h"
//#include "GeographyRateModifier.h"
//#include "RateGenerator.h"
//#include "RateMap_Biogeography.h"
//#include "RbVector.h"
//#include "Tree.h"
//#include "TypedDagNode.h"
//#include "TypedFunction.h"
//
//#include <vector>
//
//namespace RevBayesCore {
//    
//    class BiogeographyRateMapFunction : public TypedFunction<RateMap> {
//        
//    public:
//        BiogeographyRateMapFunction(size_t nc, bool fe=true, unsigned mrs=0);                                                                                             // pass in geography object??
//        virtual                                            ~BiogeographyRateMapFunction(void);                                              //!< Virtual destructor
//        
//        // public member functions
//        BiogeographyRateMapFunction*                        clone(void) const;                                                              //!< Create an independent clone
//        
//        // set parameters
//        void                                                setClockRate(const TypedDagNode< double > *r);
//        void                                                setClockRate(const TypedDagNode< RbVector< double > > *r);
//        void                                                setGeographyRateModifier(const TypedDagNode<GeographyRateModifier>* drm);
////        void                                                setCountRateModifier(const TypedDagNode<double>* rsrm);
//        void                                                setGainLossRates(const TypedDagNode< RbVector<double> >* glr);
//        void                                                setRateMatrix(const TypedDagNode<RateGenerator>* rm);
//        void                                                setRootFrequencies(const TypedDagNode< Simplex > *f);
//        void                                                update(void);
//        
//    protected:
//        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Implementation of swaping parameters
//        
//    private:
//        
//        // members
//        const TypedDagNode< double >*                       homogeneous_clock_rate;
//        const TypedDagNode< RbVector< double > >*           heterogeneous_clock_rates;
//        const TypedDagNode<RateGenerator>*                  homogeneous_rate_matrix;
//        const TypedDagNode<RbVector<RateGenerator> >*       heterogeneous_rate_matrices;
//
//        const TypedDagNode<GeographyRateModifier>*          geographyRateModifier;
////      const TypedDagNode<CountRateModifier>*              countRateModifier;
//        const TypedDagNode<Tree>*                           tau;
//        const TypedDagNode<Simplex>*                        root_frequencies;
//        
//        // geography epochs
//        // rate epochs
//        // branchwise models, etc
//        
//        bool                                                branch_heterogeneous_clock_rates;
//        bool                                                branchHeterogeneousRateMatrices;
//        bool                                                useGeographicDistance;
//        
//    };
//    
//}
//
//#endif /* defined(__rb_mlandis__BiogeographyRateMapFunction__) */
