//
//  StitchTreeFunction.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/2/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#ifndef StitchTreeFunction_h
#define StitchTreeFunction_h

#include <map>
#include <set>
#include <cstddef>
#include <vector>

#include "Tree.h"
#include "TypedFunction.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"

namespace RevBayesCore {
class DagNode;
class TopologyNode;
template <class valueType> class TypedDagNode;
    
    class StitchTreeFunction : public TypedFunction<Tree> {
        
    public:
        StitchTreeFunction(const TypedDagNode<Tree> *bt, const TypedDagNode<RbVector<Tree> >* pc, RbVector<RbVector<Taxon> > pt);
        virtual                                            ~StitchTreeFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        StitchTreeFunction*                                 clone(void) const;                                                                  //!< Create an independent clone
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        void                                                updateStitchTree(void);
        void                                                initTaxonGroups(void);
        void                                                recursivelyStitchPatchClades(TopologyNode* node, size_t& index);
        void                                                recursivelyCleanPatchClade(TopologyNode* node, TopologyNode*& newRoot, std::set<Taxon>& remainingTaxa, size_t& index, size_t patchIndex);
        void                                                recursivelyIndexPatchClade(TopologyNode* patchRoot, size_t& index, size_t patchIndex);
//        void                                                stitchPatchClade(TopologyNode* node, const Tree& patchClade);
        
        // parameters
        const TypedDagNode<Tree>*                           backboneTree;
        const TypedDagNode<RbVector<Tree> >*                patchClades;
        RbVector<RbVector<Taxon> >                          patchTaxa;
        RbVector<Taxon>                                     stitchTaxon;
        std::vector<std::set<Taxon> >                       prunedTaxa;
        size_t                                              numPatches;
        std::map<size_t, std::vector<size_t> >              stitchTreeIndex;
        std::map<Taxon, size_t>                             stitchTipIndex;
        bool                                                haveIndex;
        size_t                                              num_nodes;
        size_t                                              numTips;
        
    };
    
}

#endif /* StitchTreeFunction_h */
