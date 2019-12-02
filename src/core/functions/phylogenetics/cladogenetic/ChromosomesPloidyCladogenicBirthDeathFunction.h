//
//  ChromosomesPloidyCladogenicBirthDeathFunction.h
//
//  Created by Will Freyman on 5/22/17.
//

#ifndef ChromosomesPloidyCladogenicBirthDeathFunction__
#define ChromosomesPloidyCladogenicBirthDeathFunction__

#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <sstream>
#include <vector>
#include <map>

namespace RevBayesCore {
    
    class ChromosomesPloidyCladogenicBirthDeathFunction : public AbstractCladogenicStateFunction, public TypedFunction<CladogeneticSpeciationRateMatrix> {
        
    public:
        
        ChromosomesPloidyCladogenicBirthDeathFunction( const TypedDagNode< RbVector< double > >* sr, unsigned mc );
        virtual                                            ~ChromosomesPloidyCladogenicBirthDeathFunction(void);
        
        const static unsigned NO_CHANGE_D                           = 0;         // N -> N Diploid lineages
        const static unsigned FISSION_D                             = 1;         // N -> N + 1
        const static unsigned FUSION_D                              = 2;         // N -> N - 1
        const static unsigned POLYPLOIDIZATION_D                    = 3;         // N -> N * 2
        const static unsigned DEMIPOLYPLOIDIZATION_D                = 4;         // N -> N * 1.5
        const static unsigned NO_CHANGE_P                           = 5;         // N -> N Polyploid lineages
        const static unsigned FISSION_P                             = 6;         // N -> N + 1
        const static unsigned FUSION_P                              = 7;         // N -> N - 1
        const static unsigned POLYPLOIDIZATION_P                    = 8;         // N -> N * 2
        const static unsigned DEMIPOLYPLOIDIZATION_P                = 9;         // N -> N * 1.5
        
        // public member functions
        ChromosomesPloidyCladogenicBirthDeathFunction*              clone(void) const;
        std::map< std::vector<unsigned>, double >                   getEventMap(double t=0.0);
        const std::map< std::vector<unsigned>, double >&            getEventMap(double t=0.0) const;
        void                                                        update(void);
        
    protected:
        
        void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);
        
    private:
        
        void                                                        buildEventMap(void);
        void                                                        updateSpeciationRates(void);
        
        // members
        const TypedDagNode< RbVector<double> >*                     speciation_rates;
        unsigned                                                    max_chromo;
        unsigned                                                    num_event_types;
        std::map< std::vector<unsigned>, std::vector<unsigned> >    event_map_types;
        std::map< std::vector<unsigned>, double >                   event_map;
        std::vector< std::vector<unsigned> >                        event_map_counts;
        
    };
    
}

#endif
