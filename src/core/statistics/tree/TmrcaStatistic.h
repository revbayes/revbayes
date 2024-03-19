/**
 * @file
 * This file contains the declaration of the time of the MRCA statistic, which simply computes
 * the time of the most recent common ancestor of a clade.
 *
 * @brief Declaration of the TmrcaStatistic
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-04-13, version 1.0
 *
 * $Id: TimeTree.h 1651 2012-07-05 14:47:08Z hoehna $
 */

#ifndef TmrcaStatistic_H
#define TmrcaStatistic_H

#include <cstddef>

//#include "Statistic.h"
#include "Clade.h"
#include "TreeChangeEventListener.h"
#include "TypedFunction.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class TmrcaStatistic : public TypedFunction<double>, TreeChangeEventListener {
        
    public:
        TmrcaStatistic(const TypedDagNode<Tree> *t, const Clade &c, const bool s);                                                                                   //!< Default constructor
        virtual                                    ~TmrcaStatistic(void);                                                                  //!< Destructor
                
        // Basic utility functions
        TmrcaStatistic*                             clone(void) const;                                                                          //!< Clone object
        void                                        fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which p
        void                                        update(void);                                                                               //!< Clone the function
        
    protected:
        void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                            //!< Implementation of swaping parameters
        
    private:
        void                                        initialize(void);
        
        // members
        const TypedDagNode<Tree>*                   tree;
        Clade                                       clade;
        bool                                        initialized;
        bool                                        stem_age;
        int                                         index;
        size_t                                      taxa_count;
    };
    
}

#endif

