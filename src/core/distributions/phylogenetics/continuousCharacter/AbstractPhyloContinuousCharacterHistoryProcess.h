#ifndef AbstractPhyloContinuousCharacterHistoryProcess_H
#define AbstractPhyloContinuousCharacterHistoryProcess_H

#include <cstddef>
#include <vector>

#include "AbstractPhyloContinuousCharacterProcess.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-01-23, version 1.0
     */
    class AbstractPhyloContinuousCharacterHistoryProcess : public AbstractPhyloContinuousCharacterProcess {

    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        AbstractPhyloContinuousCharacterHistoryProcess(const TypedDagNode<Tree> *t, size_t nSites );
        virtual                                                     ~AbstractPhyloContinuousCharacterHistoryProcess(void);                                                     //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual AbstractPhyloContinuousCharacterHistoryProcess*     clone(void) const = 0;                                                                  //!< Create an independent clone
        virtual double                                              computeLnProbability(void) = 0;
        
        void                                                        redrawValue(void);
        void                                                        setRootState(const TypedDagNode< double >* s);
        void                                                        setRootState(const TypedDagNode< RbVector< double > >* s);

        
    protected:
        
        std::vector<double>                                         simulateRootCharacters(size_t n);
        virtual void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                         //!< Swap a parameter

        // members
        const TypedDagNode< double >*                               homogeneous_root_state;
        const TypedDagNode< RbVector< double > >*                   heterogeneous_root_state;

    private:
        
        double                                                      computeRootState(size_t site_index);
        
    };
    
}


#endif
