////
////  GeneralRateMapFunction.h
////  rb_mlandis
////
////  Created by Michael Landis on 4/3/14.
////  Copyright (c) 2014 Michael Landis. All rights reserved.
////
//
//#ifndef __rb_mlandis__GeneralRateMapFunction__
//#define __rb_mlandis__GeneralRateMapFunction__
//
//#include "GeographyRateModifier.h"
//#include "RateMap.h"
//#include "RateMatrix.h"
//#include "RbVector.h"
//#include "Tree.h"
//#include "TypedDagNode.h"
//#include "TypedFunction.h"
//
//#include <vector>
//
//namespace RevBayesCore {
//    
//    class GeneralRateMapFunction : public TypedFunction<RateMap> {
//        
//    public:
//        GeneralRateMapFunction(size_t ns, size_t nc);
//        GeneralRateMapFunction(const GeneralRateMapFunction &n);                                                                        //!< Copy constructor
//        virtual                                            ~GeneralRateMapFunction(void);                                               //!< Virtual destructor
//        
//        // public member functions
//        GeneralRateMapFunction*                             clone(void) const;                                                          //!< Create an independent clone
//        
//        // set parameters
//        void                                                setClockRate(const TypedDagNode< double > *r);
//        void                                                setClockRate(const TypedDagNode< RbVector< double > > *r);
//        void                                                setRateMatrix(const TypedDagNode<RateGenerator>* rm);
////        void                                                setRateMatrix(const TypedDagNode<RbVector<RateGenerator> > >* rm);
//        void                                                setRootFrequencies(const TypedDagNode< Simplex > *f);
//        void                                                update(void);
//        
//    protected:
//        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Implementation of swapping parameters
//        
//    private:
//        
//        // members
//        const TypedDagNode< double >*                               homogeneous_clock_rate;
//        const TypedDagNode< RevBayesCore::RbVector< double > >*     heterogeneous_clock_rates;
//        const TypedDagNode<RateGenerator>*                             homogeneous_rate_matrix;
//        const TypedDagNode<RbVector<RateGenerator> >*                  heterogeneous_rate_matrices;
////        const TypedDagNode<TimeTree>*                               tau;
//        const TypedDagNode<Simplex>*                        root_frequencies;
//        
//        // geography epochs
//        bool                                                branch_heterogeneous_clock_rates;
//        bool                                                branchHeterogeneousRateMatrices;
//        
//    };
//    
////    std::ostream& operator<<(std::ostream& o, const std::vector<std::vector<double> >& x);
//    
//}
//
//#endif /* defined(__rb_mlandis__GeneralRateMapFunction__) */
