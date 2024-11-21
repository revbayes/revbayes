//
//  CladogeneticProbabilityMatrixFunction.h
//
//  Created by Will Freyman on 8/1/17.
//

#ifndef CladogeneticProbabilityMatrixFunction__
#define CladogeneticProbabilityMatrixFunction__

#include <vector>
#include <cstdint>
#include <map>

#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "TypedFunction.h"

namespace RevBayesCore {
    class DagNode;
    template <class valueType> class RbVector;
    template <class valueType> class TypedDagNode;
    
    class CladogeneticProbabilityMatrixFunction : public AbstractCladogenicStateFunction, public TypedFunction<CladogeneticProbabilityMatrix> {
        
    public:
        
        CladogeneticProbabilityMatrixFunction( const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > >* events, const TypedDagNode<RevBayesCore::RbVector<double> >* probs, int n_states );
        virtual                                                                         ~CladogeneticProbabilityMatrixFunction(void);
        
        // public member functions
        CladogeneticProbabilityMatrixFunction*                                          clone(void) const;
        std::map< std::vector<unsigned>, double >                                       getEventMap(double t=0.0);
        const std::map< std::vector<unsigned>, double >&                                getEventMap(double t=0.0) const;
        void                                                                            update(void);
        
    protected:
        
        void                                                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        
        //        void                                                                            buildEventMap(void);
        void                                                                            updateProbabilitys(void);
        
        // members
        const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > >*    cladogenetic_events;
        std::map< std::vector<unsigned>, double >                                       event_map;
        unsigned                                                                        num_states;
        const TypedDagNode< RbVector<double> >*                                         probabilities;
        
    };
    
}

#endif
