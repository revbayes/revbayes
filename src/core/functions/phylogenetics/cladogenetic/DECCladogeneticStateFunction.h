#ifndef DECCladogeneticStateFunction_H
#define DECCladogeneticStateFunction_H

#include <cstddef>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <string>

#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class BranchHistory;
class DagNode;
class Simplex;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class DECCladogeneticStateFunction : public AbstractCladogenicStateFunction, public TypedFunction<CladogeneticProbabilityMatrix> {
        
    public:
        
        DECCladogeneticStateFunction(const TypedDagNode< Simplex >* ep,
                                     const TypedDagNode<RbVector<RbVector<double> > >* cg,
                                     const TypedDagNode<RbVector<RbVector<double> > >* vg,
                                     unsigned nc,
                                     unsigned ns,
                                     std::vector<std::string> et,
                                     bool epawa=true,
                                     bool wa=false,
                                     bool uv=false,
                                     unsigned mrs=0);
        
        virtual                                                 ~DECCladogeneticStateFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        DECCladogeneticStateFunction*                           clone(void) const;                                                              //!< Create an independent clone
        virtual double computeDataAugmentedCladogeneticLnProbability( const std::vector<BranchHistory*>& histories,
                                                                     size_t node_index,
                                                                     size_t left_index,
                                                                     size_t right_index ) const;
        virtual std::string simulateDataAugmentedCladogeneticState(std::vector<BranchHistory*>& histories,
                                                            size_t node_index, size_t left_index, size_t right_index) const;

        std::map< std::vector<unsigned>, double >               getEventMap(double t=0.0);
        const std::map< std::vector<unsigned>, double >&        getEventMap(double t=0.0) const;
        const std::vector<std::string>&                         getEventTypes(void) const;
        void                                                    update(void);
        
        const static unsigned MAX_NUM_AREAS = 20;
        
    protected:
        
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        void                                                    buildRanges(std::set<unsigned>& range_set, const TypedDagNode< RbVector<RbVector<double> > >* g, bool all=true);
        void                                                    buildRangesRecursively(std::set<unsigned> s, std::set<std::set<unsigned> >& r, size_t k, const TypedDagNode< RbVector<RbVector<double> > >* g, bool all=true);
        void                                                    buildBits(void);
        void                                                    buildEventMap(void);
        unsigned                                                bitsToState( const std::vector<unsigned>& b );
        std::string                                             bitsToString( const std::vector<unsigned>& b );
        std::vector<unsigned>                                   bitAllopatryComplement( const std::vector<unsigned>& mask, const std::vector<unsigned>& base );
        void                                                    bitCombinations(std::vector<std::vector<unsigned> >& comb, std::vector<unsigned> array, int i, std::vector<unsigned> accum);
        size_t                                                  computeNumStates(size_t numAreas, size_t maxRangeSize);
        unsigned                                                sumBits(const std::vector<unsigned>& b);
        void                                                    updateProbs(void);
        
        // members
        const TypedDagNode< Simplex >*                          eventProbs;
        const TypedDagNode< RbVector<RbVector<double> > >*      connectivityGraph;
        const TypedDagNode< RbVector<RbVector<double> > >*      vicarianceGraph;
        unsigned                                                numCharacters;
        unsigned                                                num_states;
        unsigned                                                numIntStates;
        unsigned                                                numRanges;
        unsigned                                                numEventTypes;
        unsigned                                                maxRangeSize;
       
        // range codes
        std::vector<std::vector<unsigned> >                     bits;
        std::map<std::vector<unsigned>, unsigned>               inverseBits;
        std::vector<std::vector<std::vector<unsigned> > >       bitsByNumOn;
        std::vector<std::vector<unsigned> >                     statesToBitsByNumOn;
        std::map< std::vector<unsigned>, unsigned>              bitsToStatesByNumOn;

        // range events: types, probs, and counts
        std::map< std::vector<unsigned>, unsigned >             eventMapTypes;
        std::map< unsigned, std::vector<unsigned> >             eventMapCounts;

        // manages simplex over event type probabilities
        std::vector<std::string>                                eventTypes;
        std::map<std::string, unsigned>                         eventStringToStateMap;
        
        // manage ranges under connectivity graph
        std::set<unsigned>                                      beforeRanges;
        std::set<unsigned>                                      afterRanges;
        
        bool                                                    eventProbsAsWeightedAverages;
        bool                                                    wideAllopatry;
        bool                                                    useVicariance;
        
    };
    
}

#endif /* defined(__revbayes_proj__DECCladogeneticStateFunction__) */
