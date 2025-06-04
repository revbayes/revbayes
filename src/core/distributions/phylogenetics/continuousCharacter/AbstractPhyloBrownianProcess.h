#ifndef AbstractPhyloBrownianProcess_H
#define AbstractPhyloBrownianProcess_H

#include <cstddef>
#include <vector>

#include "AbstractPhyloContinuousCharacterProcess.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class ContinuousTaxonData;
class Tree;
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
    class AbstractPhyloBrownianProcess : public AbstractPhyloContinuousCharacterProcess {
        
    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        AbstractPhyloBrownianProcess(const TypedDagNode<Tree> *t, size_t ns );
        virtual                                                            ~AbstractPhyloBrownianProcess(void);                                                     //!< Virtual destructor
        
        // public member functions
        // pure virtual
        virtual AbstractPhyloBrownianProcess*                               clone(void) const = 0;                                                                  //!< Create an independent clone
        virtual double                                                      computeLnProbability(void) = 0;
        
    protected:
        
        // virtual methods
        virtual void                                                        simulateRecursively(const TopologyNode& node, std::vector< ContinuousTaxonData > &t);
        virtual std::vector<double>                                         simulateRootCharacters(size_t n) = 0;
        
    };
    
}


#endif
