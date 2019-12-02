/**
 * @file
 * This file contains the declaration of the Rooted Triplet Distribution function class.
 * This class is derived from the function class and is used to
 * compute the Rooted Triplet Distribution from a set of trees.
 *
 * @brief Declaration of the Rooted Triplet Distribution function.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-07-06, version 1.0
 * @interface Function
 *
 * $Id$
 */



#ifndef RootedTripletDistributionFunction_H
#define RootedTripletDistributionFunction_H

#include "RateMatrix_GTR.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "RootedTripletDistribution.h"

#include <vector>

namespace RevBayesCore {
    
    class RootedTripletDistributionFunction : public TypedFunction<RootedTripletDistribution> {
        
    public:
        RootedTripletDistributionFunction( const TypedDagNode< RbVector< Tree > > *ts, const TypedDagNode< RbVector< std::string > > *sn );
        RootedTripletDistributionFunction( const TypedDagNode< RbVector< std::string > > *sn );
        RootedTripletDistributionFunction( const TypedDagNode< RbVector< Taxon > > *t );

        virtual                                            ~RootedTripletDistributionFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        RootedTripletDistributionFunction*                  clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
        void                                                setTrees(TypedDagNode< RbVector< Tree > >* ts); //!< Set the trees in the distribution

        void                                                setRecordBranchLengths(  TypedDagNode< bool >* b) ; //!< Set whether we record some branch lengths information
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        void                                                extractTriplets( Tree& t );                           //!< Getting the triplets from a tree and adding them to the tripletDistribution
        void                                                extractTriplets(  );             //!< Getting the triplets from a tree and adding them to the tripletDistribution
        void                                                populateTripletDistribution ( TopologyNode* node, std::vector< size_t >& allTips ) ;
        void                                                addAllTriplets(std::vector< size_t >& leftTips, std::vector< size_t >& rightTips) ; //!< Get all rooted triplets given vectors of left and right tips
        void                                                addAllTripletsOneWay( std::vector< size_t >& leftTips, std::vector< size_t >& rightTips, size_t leftSize,size_t rightSize ); //!< Get rooted triplets given vectors of left and right tips, one way only
        

        // members

        const TypedDagNode< RbVector<Tree> >*               trees;
        const TypedDagNode< RbVector< Taxon > >*            taxa;
        const TypedDagNode< RbVector< std::string > >*      species;
        const TypedDagNode< bool >*                         recordBranchLengths;
        size_t                                              numberOfTrees;
        std::set< Taxon >                                   taxaSet;
        std::set< std::string >                             speciesSet;
        std::map< Taxon, size_t >                           taxonToIndex;
        std::map< std::string, size_t >                     speciesToIndex;
        bool                                                speciesOnly;
        RootedTripletDistribution*                          rtd;

    };
    
}

#endif
