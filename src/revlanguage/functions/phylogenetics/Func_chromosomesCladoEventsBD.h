//
//  Func_chromosomesCladoEventsBD.h
//
//  Created by Will Freyman on 6/22/16.
//

#ifndef Func_chromosomesCladoEventsBD_H
#define Func_chromosomesCladoEventsBD_H


#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlTypedFunction.h"

#include <map>
#include <string>

namespace RevLanguage {
    
    class Func_chromosomesCladoEventsBD : public TypedFunction<CladogeneticSpeciationRateMatrix> {
        
    public:
        Func_chromosomesCladoEventsBD( void );
        

        Func_chromosomesCladoEventsBD*                                                  clone(void) const;               //!< Clone the object
        static const std::string&                                                       getClassType(void);              //!< Get Rev type
        static const TypeSpec&                                                          getClassTypeSpec(void);          //!< Get class type spec
        std::string                                                                     getFunctionName(void) const;     //!< Get the primary name of the function in Rev
        const TypeSpec&                                                                 getTypeSpec(void) const;         //!< Get the type spec of the instance
        

        RevBayesCore::TypedFunction< RevBayesCore::CladogeneticSpeciationRateMatrix>*   createFunction(void) const;      //!< Create internal function object
        const ArgumentRules&                                                            getArgumentRules(void) const;    //!< Get argument rules
        
    };
    
}




#endif
