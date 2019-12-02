#ifndef Func_simCompleteTree_H
#define Func_simCompleteTree_H

#include "Procedure.h"
#include "Tree.h"
#include "TopologyNode.h"

namespace RevLanguage {
    
    /**
     * @brief Rev function to simumlate a balanced or catapilar tree.
     *
     * This procedure simulates a simple tree.
     * The tree is either totally balanced, or a caterpillar tree.
     * The internode distance is fix to 1.0.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-02-27
     *
     */
    class Func_simCompleteTree : public Procedure {
        
    public:
        Func_simCompleteTree( void );
        
        // Basic utility functions
        Func_simCompleteTree*                       clone(void) const;                                          //!< Clone object
        static const std::string&                   getClassType(void);                                         //!< Get Rev simTree
        static const TypeSpec&                      getClassTypeSpec(void);                                     //!< Get class simTree spec
        std::string                                 getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                             getTypeSpec(void) const;                                    //!< Get language simTree of the object
        
        // Func_source functions
        const ArgumentRules&                        getArgumentRules(void) const;                               //!< Get argument rules
        const TypeSpec&                             getReturnType(void) const;                                  //!< Get simTree of return val
        
        RevPtr<RevVariable>                         execute(void);                                              //!< Execute function
        
    private:
        
        void                                        buildRandomBinaryTree(std::vector<RevBayesCore::TopologyNode*> &tips, std::vector<double> &ages);
        
    };
    
}

#endif

