//
//  GeoSSECladogeneticBirthDeathFunction.hpp
//  revbayes-proj
//
//  Created by Michael R May on 12/14/22.
//  Copyright © 2022 Michael R May. All rights reserved.
//

#ifndef GeoSSECladogeneticBirthDeathFunction__
#define GeoSSECladogeneticBirthDeathFunction__

#include <stddef.h>
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
    
    class GeoSSECladogeneticBirthDeathFunction : public AbstractCladogenicStateFunction, public TypedFunction<CladogeneticSpeciationRateMatrix> {
        
    public:
        
        GeoSSECladogeneticBirthDeathFunction(const TypedDagNode<RbVector<double> >* sr, const TypedDagNode< RbVector<double> >* ar);
        virtual                                                     ~GeoSSECladogeneticBirthDeathFunction(void);
        
        const static unsigned NUM_CLADO_EVENT_TYPES                 = 3;
        const static unsigned SYMPATRY                              = 0;         // A  -> A or AB -> AB|A
        const static unsigned ALLOPATRY                             = 1;         // AB -> A|B
        const static unsigned JUMP_DISPERSAL                        = 2;         // A  -> A|B
        
        // public member functions
        GeoSSECladogeneticBirthDeathFunction*                       clone(void) const;
        
        std::map< std::vector<unsigned>, double >                   getEventMap(double t=0.0);
        const std::map< std::vector<unsigned>, double >&            getEventMap(double t=0.0) const;        
        void                                                        setJumpRates(const TypedDagNode< RbVector< RbVector<double> > >* jr);
        void                                                        update(void);

    protected:
        
        void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        unsigned                                                    bitsToState( const std::vector<unsigned>& b );
        std::string                                                 bitsToString( const std::vector<unsigned>& b );
        std::vector<unsigned>                                       bitAllopatryComplement( const std::vector<unsigned>& mask, const std::vector<unsigned>& base );
        void                                                        bitCombinations(std::vector<std::vector<unsigned> >& comb, std::vector<unsigned> array, int i, std::vector<unsigned> accum);
        void                                                        buildBits(void);
        void                                                        buildEventMap(void);
        void                                                        buildRanges(std::set<unsigned>& range_set, bool all=true);
        void                                                        buildRangesRecursively(std::set<unsigned> s, std::set<std::set<unsigned> >& r, bool all=true);
        size_t                                                      computeNumStates(size_t numAreas, size_t maxRangeSize);
        void                                                        printEventMap();
        unsigned                                                    sumBits(const std::vector<unsigned>& b);
        void                                                        updateEventMapCutsetWeights(void);
        
        // parameters
        const TypedDagNode< RbVector<double> >*                     sympatryRates;
        const TypedDagNode< RbVector<double> >*                     allopatryRates;
        bool                                                        hasJumps;
        const TypedDagNode<RbVector<RbVector<double> > >*           jumpRates;

        // dimensions
        unsigned                                                    numCharacters;
        unsigned                                                    numStates;
        unsigned                                                    numIntStates;
        unsigned                                                    numRanges;
        
        // model settings
        
        // range codes
        std::vector<std::vector<unsigned> >                         bits;
        std::map<std::vector<unsigned>, unsigned>                   inverseBits;
        std::vector<std::vector<std::vector<unsigned> > >           bitsByNumOn;
        std::vector<std::vector<unsigned> >                         statesToBitsByNumOn;
        std::vector<std::set<unsigned> >                            statesToBitsetsByNumOn;
        std::map< std::vector<unsigned>, unsigned>                  bitsToStatesByNumOn;
        std::set<unsigned>                                          ranges;
        
        // event maps
        std::map< std::vector<unsigned>, double >                   eventMap;
        std::map< std::vector<unsigned>, unsigned>                  eventMapArea;
        std::map< std::vector<unsigned>, unsigned>                  eventMapTypes;
        std::map< unsigned, std::vector<unsigned> >                 eventMapCounts;
        
        // MJL: eventually, deprecate this stuff
        // manages string-based simplex mapping??
        std::vector<std::string>                                    eventTypes;
        std::map<std::string, unsigned>                             eventStringToStateMap;

    };
    
}

#endif
