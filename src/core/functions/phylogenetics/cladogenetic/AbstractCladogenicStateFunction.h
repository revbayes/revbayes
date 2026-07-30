#ifndef AbstractCladogenicStateFunction_h
#define AbstractCladogenicStateFunction_h

#include "BranchHistory.h"
#include "CharacterHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "TypedDagNode.h"

#include <vector>

namespace RevBayesCore {
    
    class AbstractCladogenicStateFunction {
    
    public:
        
        virtual std::map< std::vector<unsigned>, double >                       getEventMap(double t=0.0) = 0;
        virtual const std::map< std::vector<unsigned>, double >&                getEventMap(double t=0.0) const = 0;
        virtual double computeDataAugmentedCladogeneticLnProbability( const CharacterHistoryDiscrete& histories,
                                                                      size_t node_index,
                                                                      size_t left_index,
                                                                      size_t right_index ) const { return 0.0; };
        virtual std::string simulateDataAugmentedCladogeneticState(CharacterHistoryDiscrete& histories,
                                                            size_t node_index, size_t left_index, size_t right_index) const { return ""; };
        
	virtual ~AbstractCladogenicStateFunction() {};

    protected:
        
        virtual void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP) = 0;
        
    };
    
}

#endif

