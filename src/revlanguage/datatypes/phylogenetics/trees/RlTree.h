#ifndef RlTree_H
#define RlTree_H

#include "ModelObject.h"
#include "Tree.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    
    /**
     * @file
     * This file contains the declaration of a Tree.
     *
     * @brief Declaration of Tree
     *
     * (c) Copyright 2009-
     * @date Last modified: $Date: $
     * @author The RevBayes Development Core Team
     * @license GPL version 3
     *
     *
     */
    class Tree : public ModelObject<RevBayesCore::Tree> {
        
    public:
                                            Tree(void);                                                                         //!< Constructor requires character type
                                            Tree(RevBayesCore::Tree *v);                                                        //!< Constructor requires character type
                                            Tree(const RevBayesCore::Tree &v);                                                  //!< Constructor requires character type
                                            Tree(RevBayesCore::TypedDagNode<RevBayesCore::Tree> *n);                            //!< Constructor requires character type
        typedef RevBayesCore::Tree          valueType;
        
        // (pure) virtual methods
        virtual Tree*                       clone(void) const;                                                                  //!< Clone object
        
        // Basic utility functions
        static const std::string&           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                             //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                            //!< Get language type of the object
        
        // Member method functions
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f); //!< Map member methods to internal functions
        virtual std::string                 getGuiName(void) { return "Topology"; }
        virtual std::string                 getGuiUnicodeSymbol(void) { return "\\tau"; }
        std::string                         getGuiInfo(void) { return ""; }
        
    protected:
        void                                initMethods(void);
    };
    
}

#endif
