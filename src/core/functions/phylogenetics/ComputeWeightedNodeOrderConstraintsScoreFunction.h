//
//  computeWeightedNodeOrderConstraintsScoreFunction.h
//  revbayes-proj
//
//  Created by Bastien Boussau on 16/11/15.
//  Copyright (c) 2015 Bastien Boussau. All rights reserved.
//

#ifndef computeWeightedNodeOrderConstraintsScoreFunction_H
#define computeWeightedNodeOrderConstraintsScoreFunction_H

#include "RelativeNodeAgeWeightedConstraints.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <set>

#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/unordered_set.hpp>

namespace RevBayesCore {
    
    class computeWeightedNodeOrderConstraintsScoreFunction : public TypedFunction<double> {
        
    public:
        computeWeightedNodeOrderConstraintsScoreFunction(const TypedDagNode<Tree> *t, const TypedDagNode<RelativeNodeAgeWeightedConstraints> *co, TypedDagNode< double>* bet);
        virtual                                            ~computeWeightedNodeOrderConstraintsScoreFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        computeWeightedNodeOrderConstraintsScoreFunction*                                  clone(void) const;                                                                  //!< Create an independent clone
        void                                                keep(DagNode* affecter);
        void                                                restore(DagNode *restorer);
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                touch(DagNode *toucher );
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Tree>*                           tau;
        const TypedDagNode<RelativeNodeAgeWeightedConstraints>*         constraints;
        const TypedDagNode< double>* beta;
        std::map<std::pair<std::string, std::string>, double >  nodeAges;
// @Bastien: We currently do not support c++11 because of compile issues on some older machines.
//        std::unordered_set< std::pair < std::string, std::string >, boost::hash<std::pair<std::string, std::string> > > constrainedNodes;
        std::set< std::pair < std::string, std::string > > constrainedNodes;

        // functions
        void updateSetOfConstrainedNodes();
        void updateMapOfNodeAges();
        double computeWeightedNodeAgeConstraintsScore();

    };
    
}

#endif
