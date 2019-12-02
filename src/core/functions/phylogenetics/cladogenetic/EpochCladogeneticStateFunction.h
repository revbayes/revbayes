//
//  EpochCladogeneticStateFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef EpochCladogeneticStateFunction_h
#define EpochCladogeneticStateFunction_h

#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticProbabilityMatrix_Epoch.h"
#include "RbVector.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <sstream>
#include <vector>
#include <map>

namespace RevBayesCore {
    
    class EpochCladogeneticStateFunction : public AbstractCladogenicStateFunction, public TypedFunction<CladogeneticProbabilityMatrix> {
        
    public:
        
        EpochCladogeneticStateFunction( const TypedDagNode< RbVector< double > >* et, const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* cp, unsigned nc, unsigned ns);
        virtual                                                 ~EpochCladogeneticStateFunction(void);                                              //!< Virtual destructor
        
        // public member functions
        EpochCladogeneticStateFunction*                         clone(void) const;                     //!< Create an independent clone
        std::map< std::vector<unsigned>, double >               getEventMap(double t=0.0);
        const std::map< std::vector<unsigned>, double >&        getEventMap(double t=0.0) const;
        
        void                                                    update(void);
        
    protected:
        
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
//        size_t                                                  findEpochIndex( double t ) const;
//        void                                                    buildEventMap(void);
        
        // members
        const TypedDagNode< RbVector<double> >*                 epochTimes;
        const TypedDagNode< RbVector<CladogeneticProbabilityMatrix> >*             cladoProbs;
        unsigned                                                numCharacters;
        unsigned                                                numStates;
        unsigned                                                numEpochs;

    };
    
}

#endif /* EpochCladogeneticStateFunction_h */
