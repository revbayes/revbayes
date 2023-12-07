//
//  ChromosomesCladogenicBirthDeathFunction.h
//
//  Created by Will Freyman on 6/22/16.
//

#ifndef ChromosomesCladogenicBirthDeathFunction__
#define ChromosomesCladogenicBirthDeathFunction__

#include <cstddef>
#include <vector>
#include <map>

#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class BranchHistory;
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class ChromosomesCladogenicBirthDeathFunction : public AbstractCladogenicStateFunction, public TypedFunction<CladogeneticSpeciationRateMatrix> {
        
    public:
        
        ChromosomesCladogenicBirthDeathFunction( const TypedDagNode< RbVector< double > >* sr, unsigned mc );
        virtual                                                     ~ChromosomesCladogenicBirthDeathFunction(void);
        
        const static unsigned NO_CHANGE                             = 0;         // N -> N
        const static unsigned FISSION                               = 1;         // N -> N + 1
        const static unsigned FUSION                                = 2;         // N -> N - 1
        const static unsigned POLYPLOIDIZATION                      = 3;         // N -> N * 2
        const static unsigned DEMIPOLYPLOIDIZATION                  = 4;         // N -> N * 1.5
        
        // public member functions
        virtual double computeDataAugmentedCladogeneticLnProbability( const std::vector<BranchHistory*>& histories,
                                                                     size_t node_index,
                                                                     size_t left_index,
                                                                     size_t right_index ) const;

        ChromosomesCladogenicBirthDeathFunction*                    clone(void) const;
        std::map< std::vector<unsigned>, double >                   getEventMap(double t=0.0);
        const std::map< std::vector<unsigned>, double >&            getEventMap(double t=0.0) const;
        void                                                        setRateMultipliers(const TypedDagNode< RbVector< double > >* rm);
        void                                                        update(void);
        
    protected:
        
        void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        
        void                                                        buildEventMap(void);
        void                                                        updateSpeciationRates(void);
        
        // members
        const TypedDagNode< RbVector<double> >*                     speciationRates;
        const TypedDagNode< RbVector<double> >*                     hiddenRateMultipliers;
        unsigned                                                    maxChromo;
        unsigned                                                    numEventTypes;
        std::map< std::vector<unsigned>, std::vector<unsigned> >    eventMapTypes;
        std::map< std::vector<unsigned>, double >                   eventMap;
        std::vector< std::vector<unsigned> >                        eventMapCounts;
        bool                                                        use_hidden_rate;
        
    };
    
}

#endif
