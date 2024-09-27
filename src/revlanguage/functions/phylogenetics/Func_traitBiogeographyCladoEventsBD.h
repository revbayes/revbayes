//
//  Func_traitBiogeographyCladoEventsBD.h
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/27/24.
//

#ifndef Func_traitBiogeographyCladoEventsBD_h
#define Func_traitBiogeographyCladoEventsBD_h

#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_traitBiogeographyCladoEventsBD : public TypedFunction<CladogeneticSpeciationRateMatrix> {
        
    public:
        Func_traitBiogeographyCladoEventsBD( void );
        
        
        Func_traitBiogeographyCladoEventsBD*                                            clone(void) const;               //!< Clone the object
        static const std::string&                                                       getClassType(void);              //!< Get Rev type
        static const TypeSpec&                                                          getClassTypeSpec(void);          //!< Get class type spec
        std::string                                                                     getFunctionName(void) const;     //!< Get the primary name of the function in Rev
        const TypeSpec&                                                                 getTypeSpec(void) const;         //!< Get the type spec of the instance
        
        
        RevBayesCore::TypedFunction< RevBayesCore::CladogeneticSpeciationRateMatrix>*   createFunction(void) const;      //!< Create internal function object
        const ArgumentRules&                                                            getArgumentRules(void) const;    //!< Get argument rules
        
    };
    
}

#endif /* Func_traitBiogeographyCladoEventsBD_hpp */
