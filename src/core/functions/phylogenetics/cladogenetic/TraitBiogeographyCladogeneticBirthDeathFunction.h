//
//  TraitBiogeographyCladogeneticBirthDeathFunction.hpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/27/24.
//

#ifndef TraitBiogeographyCladogeneticBirthDeathFunction_h
#define TraitBiogeographyCladogeneticBirthDeathFunction_h

#include <cstddef>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>

#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "CladogeneticProbabilityMatrix.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class BranchHistory;
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class TraitBiogeographyCladogeneticBirthDeathFunction :
        public AbstractCladogenicStateFunction,
        public TypedFunction<CladogeneticSpeciationRateMatrix> {
    
        // type definitions
        typedef std::vector<unsigned> TraitBits;
        typedef std::vector<unsigned> RangeBits;
        typedef std::set<unsigned> RangeBitset;
        typedef std::set<unsigned> TraitBitset;
        typedef std::vector<unsigned> StateTriplet;
//        typedef std::vector<unsigned> TraitDoublet;
        typedef std::pair<RangeBits, TraitBits> CompositeBits;
        typedef std::pair<RangeBitset, TraitBitset> CompositeBitset;
            
        public:
            
            TraitBiogeographyCladogeneticBirthDeathFunction(const TypedDagNode<RbVector<RbVector<double> > >* rw,
                                                            const TypedDagNode<RbVector<RbVector<double> > >* rb,
                                                            const TypedDagNode<RbVector<RbVector<RbVector<double> > > >* mw,
                                                            const TypedDagNode<RbVector<RbVector<RbVector<RbVector<double> > > > >* mb,
                                                            unsigned mrs, unsigned msss, bool nss, std::string ct="cutset");

            virtual                                                     ~TraitBiogeographyCladogeneticBirthDeathFunction(void);
            
            
            // public member functions
            virtual double computeDataAugmentedCladogeneticLnProbability(const std::vector<BranchHistory*>& histories, size_t node_index, size_t left_index, size_t right_index ) const;
            
            TraitBiogeographyCladogeneticBirthDeathFunction*  clone(void) const;
            
            std::map<StateTriplet, double >                   getEventMap(double t=0.0);
            const std::map<StateTriplet, double >&            getEventMap(double t=0.0) const;
           
            void                                              update(void);
            
            const static unsigned NUM_CLADO_EVENT_TYPES = 3;
            
            const static unsigned WITHIN_SPECIATION     = 0; // A  -> A or AB -> AB|A
            const static unsigned BETWEEN_SPECIATION    = 1; // AB -> A|B
            const static unsigned FOUNDER_SPECIATION    = 2; // A  -> A|B
        
        protected:
            
            void                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);
            
        private:
            unsigned                compositeBitsToState(const CompositeBits& b);
            std::string             compositeBitsToString(const CompositeBits& b);
            std::vector<unsigned>   rangeBitComplement(const RangeBits& mask,
                                                       const RangeBits& base);
            void                    rangeBitCombinations(std::vector<RangeBits>& comb,
                                                         RangeBits array, int i,
                                                         RangeBits accum);
            void                    buildStateSpace(void);
            void                    buildBuddingRegions(void);
            void                    buildCutsets(void);
            void                    buildEventMap(void);
            void                    buildEventMapFactors(void);
            void                    buildRanges(std::set<unsigned>& range_set,
                                                const RbVector<RbVector<double> >& g,
                                                bool all=true);
            void                    buildRangesRecursively(std::set<unsigned> s,
                                                           std::set<std::set<unsigned> >& r,
                                                           size_t k,
                                                           const RbVector<RbVector<double> >& g,
                                                           bool all=true);
            double                  computeCutsetScore(StateTriplet idx, unsigned et, size_t trait_idx, size_t val_idx);
//            double                  computeModularityScore(std::vector<unsigned> idx, unsigned et);
            size_t                  computeNumStates(size_t numAreas, size_t maxRangeSize);
            void                    printEventMap(std::map<StateTriplet, double > x);
            unsigned                sumRangeBits(const RangeBits& b);
            void                    updateEventMapCutsetWeights(void);
            
            // parameters
            const TypedDagNode<RbVector<RbVector<double> > >*                         rho_w;  // within-region speciation base rate
            const TypedDagNode<RbVector<RbVector<double> > >*                         rho_b;  // between-region speciation base rate
            const TypedDagNode<RbVector<RbVector<RbVector<double> > > >*              m_w;    // within-region speciation relative rates
                                                                                              // dim1: trait-chars, dim2:                 trait-vals, dim3: regions
            const TypedDagNode<RbVector<RbVector<RbVector<RbVector<double> > > > >*   m_b;    // between-region speciation relative rates
                                                                                              // dim1: trait-chars, dim2: trait-vals,
                                                                                              // dim3: regions, dim4: regions
            
            // dimensions
            size_t        numTraits;        // number of traits
            size_t        numTraitValues;   // number of values per trait
            size_t        numTraitSets;     // number of trait-value pairs
            size_t        numRegions;       // number of discrete regions
            size_t        numIntStates;     // number of integer states
                                            // accounts for max range size, etc.
            size_t        numRanges;        // number of presence-absence ranges
            size_t        numCompositeStates;
            size_t        maxRangeSize;
            size_t        maxSubrangeSplitSize;
            
            // model settings
            bool            use_cutset_mean;
            bool            normalize_split_scores;
            std::string     connectivityType;
            
            // range codes
            std::vector<RangeBits>                         rangeBits;
            std::vector<TraitBits>                         traitBits;
            std::vector<std::vector<RangeBits> >           rangeBitsByNumOn;
            std::vector<RangeBits>                         statesToRangeBitsByNumOn;
            std::vector<RangeBitset>                       statesToRangeBitsetsByNumOn;
            std::vector<TraitBits>                         statesToTraitBitsByNumOn;
            std::vector<TraitBitset>                       statesToTraitBitsetsByNumOn;
            std::vector<CompositeBits>                     statesToCompositeBitsByNumOn;
            std::vector<CompositeBitset>                   statesToCompositeBitsetsByNumOn;
            std::map< CompositeBits, unsigned>             compositeBitsToStatesByNumOn;
            std::set<unsigned>                             ranges;
            std::set<unsigned>                             traits;
            
            // event maps
            size_t                               numEventTypes;
            std::map< StateTriplet, double >     eventMap;                // Key is state triplet, value is absolute event rate
//            std::map< TraitBits, std::map<StateTriplet, double > > eventMapByTraits;                // Key is state triplet, value is absolute event rate
            std::map< StateTriplet, unsigned>    eventMapTypes;           // Key is state triplet, value is event type (e.g. within vs. between speciation)
            std::map< size_t, std::vector<unsigned> > eventMapCounts;   // Key is integer-valued state, value is number of ways to speciate
            std::map< StateTriplet, std::vector< std::vector<unsigned> > > eventMapCutsets; // returns the vector of cut edges for a given left/right split
            std::map< StateTriplet, unsigned >   eventMapBuddingRegions; // returns the vector of cut edges for a given left/right split
//            std::map< TraitDoublet, std::vector<StateTriplet> >   eventMapTraitToRangePattern;
            std::vector<std::vector< std::map< StateTriplet, double > > >    eventMapFactors;         // Key is state triplet, value is *unnormalized* relative speciation rate
            std::vector<std::vector< std::map< StateTriplet, double > > >    eventMapWeights;         // Key is state triplet, value is *normalized* relative speciation rate
            
            // MJL: eventually, deprecate this stuff
            // manages string-based simplex mapping??
            std::vector<std::string>                                    eventTypes;
            std::map<std::string, unsigned>                             eventStringToStateMap;

    };
    
}

#endif /* TraitBiogeographyCladogeneticBirthDeathFunction_h */
