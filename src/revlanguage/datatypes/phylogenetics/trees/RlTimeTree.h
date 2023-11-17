/**
 * @file
 * This file contains the declaration of a TimeTree, which is
 * class that holds a tree with nodes having times attached.
 *
 * @brief Declaration of TimeTree
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef RlTimeTree_H
#define RlTimeTree_H

#include "ModelObject.h"
#include "RlTree.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    class TimeTree : public Tree {
        
    public:
        TimeTree(void);                                                                                                             //!< Constructor requires character type
        TimeTree(RevBayesCore::Tree *v);                                                                                        //!< Constructor requires character type
        TimeTree(const RevBayesCore::Tree &v);                                                                                  //!< Constructor requires character type
        TimeTree(RevBayesCore::TypedDagNode<RevBayesCore::Tree> *n);                                                            //!< Constructor requires character type
        
        // Basic utility functions
        TimeTree*                           clone(void) const;                                                                      //!< Clone object
        static const std::string&           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object

        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);     //!< Map member methods to internal functions
        
        
    protected:
        
        void                                initMethods(void);
        void                                recursiveSetAge(RevBayesCore::Tree& tree, size_t index, double age, double diff);

    };
    
}

#endif
