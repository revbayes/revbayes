
#ifndef RlClade_H
#define RlClade_H

#include "Clade.h"
#include "ModelObject.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    /**
    * @file
    * This file contains the declaration of a Clade, which is
    * the class that contains the names of the taxa that build the clade.
    *
    * @brief Declaration of Clade
    *
    * (c) Copyright 2009-
    * @date Last modified: $Date: $
    * @author The RevBayes Development Core Team
    * @license GPL version 3
    *
    * 
    */
    class Clade : public ModelObject<RevBayesCore::Clade> {
        
        public:
                                                    Clade(void);                                                                        //!< Constructor requires character type
                                                    Clade(RevBayesCore::Clade *v);                                                      //!< Constructor requires character type
                                                    Clade(const RevBayesCore::Clade &v);                                                //!< Constructor requires character type
                                                    Clade(RevBayesCore::TypedDagNode<RevBayesCore::Clade> *n);                          //!< Constructor requires character type
            
            typedef RevBayesCore::Clade             valueType;
            
            // Basic utility functions
            Clade*                                  clone(void) const;                                                                  //!< Clone object
            void                                    constructInternalObject(void);                                                      //!< We construct the a new internal MCMC object.
            virtual RevPtr<RevVariable>             executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f); //!< Map member methods to internal functions
            static const std::string&               getClassType(void);                                                                 //!< Get Rev type
            static const TypeSpec&                  getClassTypeSpec(void);                                                             //!< Get class type spec
            std::string                             getConstructorFunctionName(void) const;                                             //!< Get the name used for the constructor function in Rev.
            const MemberRules&                      getParameterRules(void) const;                                                      //!< Get member rules (const)
            const TypeSpec&                         getTypeSpec(void) const;                                                            //!< Get language type of the object
            std::string                             getGuiName(void) { return "Clade"; }
            std::string                             getGuiUnicodeSymbol(void) { return "\\tau(S)"; }
            std::string                             getGuiInfo(void) { return ""; }
            
        protected:
            void                                    initMethods(void);
            void                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
        
            std::vector<RevPtr<const RevVariable> > names;
            std::vector<RevPtr<const RevVariable> > names_vector;
            std::vector<RevPtr<const RevVariable> > clades;
            RevPtr<const RevVariable>               age;
            RevPtr<const RevVariable>               missing;
            RevPtr<const RevVariable>               is_negative_constraint;
            RevPtr<const RevVariable>               optional_match;
            //        RevPtr<const RevVariable>               optional_constraints;

    };
    
}

#endif
