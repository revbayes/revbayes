/**
 * @file
 * This file contains the declaration and implementation of the PhyloDiversity function
 * which calculates the sum of all branches on the tree needed to connect the sampled taxa.
 *
 * @brief Declaration and implementation of PhyloDiversity
 *
 * (c) Copyright 2009- under GPL version 3
 * @author Will Freyman
 * @license GPL version 3
 * @version 1.0
 *
 */

#ifndef PhyloDiversityFunction_H
#define PhyloDiversityFunction_H


#include <cstddef>
#include <vector>

#include "Clade.h"
#include "TypedFunction.h"


namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class PhyloDiversityFunction : public TypedFunction<double> {
        
    public:
        PhyloDiversityFunction(const TypedDagNode<Tree> *t, const Clade &c, const bool i, const TypedDagNode< RbVector< double > > *w);
        virtual                                            ~PhyloDiversityFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        PhyloDiversityFunction*                             clone(void) const;                                                                  //!< Create an independent clone
        void                                                update(void);
        
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        double                                              calculateBranchWeights( size_t j );
        double                                              sumPDforNode(size_t j, std::vector<size_t>* nodesVisited, size_t stopIndex);
        
        // members
        std::vector<double>                                 branchWeights;
        bool                                                includeRoot;
        size_t                                              num_taxa;
        const Clade                                         sample;
        const TypedDagNode<Tree>*                           tau;
        const TypedDagNode< RbVector< double > >*           tipWeights;
        bool                                                weightedPD;
    };
    
}

#endif
