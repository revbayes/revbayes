/**
 * @file
 * This file contains the declaration of the Maximum Tree function class.
 * This class is derived from the function class and is used to
 * compute the maximum species tree (Liu, 2006) from a set of gene trees.
 * This tree is consistent under the multispecies coalescent with a single effective population size.
 *
 * @brief Declaration of the Maximum Tree function.
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



#ifndef MaximumTreeFunction_H
#define MaximumTreeFunction_H

#include <vector>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>
#include <utility>

#include "Tree.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
    
    class MaximumTreeFunction : public TypedFunction<Tree> {
        
    public:
        MaximumTreeFunction( const TypedDagNode< RbVector<Tree> > *ts );
        virtual                                            ~MaximumTreeFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        MaximumTreeFunction*                                clone(void) const;                                                              //!< Create an independent clone
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        void                                                getMinDepthMatrix (  );
        Tree*                                               getSpeciesTreeFromMinDepths (  ) ;

        // members
        
        const TypedDagNode< RbVector<Tree> >*               trees;
        size_t                                              numSpecies;
        std::vector<double>                                 depthMatrix;
        std::vector< std::pair<std::string, std::string> >  speciesPairSets;
        std::vector<std::string>                            speciesNamesV;
        std::map < std::string, size_t >                    speciesToIndex;

    };
    
}

#endif
