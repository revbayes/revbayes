//
//  CladogeneticProbabilityMatrix_Epoch_Epoch.h
//  revbayes-proj
//
//  Created by Michael Landis on 10/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef CladogeneticProbabilityMatrix_Epoch_Epoch_h
#define CladogeneticProbabilityMatrix_Epoch_Epoch_h

#include <cstddef>
#include <iosfwd>
#include <map>
#include <vector>

#include "CladogeneticProbabilityMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore {
class Assignable;
    
    //    class TransitionProbabilityMatrix;
    
    class CladogeneticProbabilityMatrix_Epoch : public CladogeneticProbabilityMatrix {
        
    public:
        CladogeneticProbabilityMatrix_Epoch(void);          //!< Construct rate matrix with
        CladogeneticProbabilityMatrix_Epoch(size_t n);      //!< Construct rate matrix with n states
        virtual                                             ~CladogeneticProbabilityMatrix_Epoch(void);
        
        // pure virtual methods
        virtual CladogeneticProbabilityMatrix_Epoch*        clone(void) const;
        virtual void                                        initFromString( const std::string &s );
        
        // virtual methods that may need to overwritten
        virtual void                                        update(void);
        std::map<std::vector<unsigned>, double>             getEventMap(double t=0.0);
        const std::map<std::vector<unsigned>, double>&      getEventMap(double t=0.0) const;
        const RbVector<double>&                             getEpochTimes(void) const;                                                                //!< Return the epoch times
        const RbVector<CladogeneticProbabilityMatrix>&      getCladogeneticProbabilityMatrix(void) const;
        const CladogeneticProbabilityMatrix&                getCladogeneticProbabilityMatrix(double t) const;
        void                                                setEventMap(const std::map<std::vector<unsigned>,double>& m, double t);
        void                                                setEpochCladogeneticProbabilityMatrix(const RbVector<CladogeneticProbabilityMatrix>& cp); //!< Update all epoch matrices
        void                                                setEpochCladogeneticProbabilityMatrix(const CladogeneticProbabilityMatrix& cp, double t); //!< Update only epoch matrix at time t
        void                                                setEpochTimes(const RbVector<double> &t);                                                 //!< Directly set the epoch times
        
	json                                                toJSON() const;
        virtual void                                        printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const;            //!< print object for user (in user-formatted way)
        
    protected:
        
        size_t                                              findEpochIndex( double t ) const;
        
        // parameters
        RbVector<CladogeneticProbabilityMatrix>             epochCladogeneticProbabilityMatrices;
        RbVector<double>                                    epochTimes;
        
        // helper variables
        size_t                                              numEpochs;
        bool                                                needs_update;
        
    };
    
};
#endif /* CladogeneticProbabilityMatrix_Epoch_Epoch_h */
