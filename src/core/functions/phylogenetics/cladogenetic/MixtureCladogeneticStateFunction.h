#ifndef MixtureCladogeneticStateFunction_h
#define MixtureCladogeneticStateFunction_h


#include <cstddef>
#include <vector>
#include <map>

#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
class Simplex;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class MixtureCladogeneticStateFunction : public AbstractCladogenicStateFunction, public TypedFunction<CladogeneticProbabilityMatrix> {
        
    public:
        
        MixtureCladogeneticStateFunction( const TypedDagNode< Simplex >* mw, const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* cp, unsigned nc, unsigned ns);
        virtual                                                 ~MixtureCladogeneticStateFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        MixtureCladogeneticStateFunction*                       clone(void) const;                                                              //!< Create an independent clone
        std::map< std::vector<unsigned>, double >               getEventMap(double t=0.0);
        const std::map< std::vector<unsigned>, double >&        getEventMap(double t=0.0) const;
        void                                                    update(void);
        
    protected:
        
        void                                                    swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        void                                                    buildEventMap(void);
        
        // members
        const TypedDagNode< Simplex >*                          mixtureWeights;
        const TypedDagNode< RbVector<CladogeneticProbabilityMatrix> >*             cladoProbs;
        unsigned                                                numCharacters;
        unsigned                                                numStates;
        size_t                                                  numComponents;
        
        std::map< std::vector<unsigned>, unsigned >             eventMapTypes;
//        std::map< std::vector<unsigned>, double >               eventMapProbs;
        std::vector< std::vector<unsigned> >                    eventMapCounts;
        
    };
    
}


#endif /* MixtureCladogeneticStateFunction_h */
