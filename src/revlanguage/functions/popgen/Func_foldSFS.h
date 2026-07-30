#ifndef Func_foldSFS_H
#define Func_foldSFS_H

#include "ModelVector.h"
#include "Natural.h"
#include "RlTypedFunction.h"

namespace RevLanguage {


    /**
     * The RevLanguage wrapper of the SFS folding function.
     *
     * Folds an unfolded site frequency spectrum (SFS) by combining
     * complementary allele-frequency classes.
     */
    class Func_foldSFS : public TypedFunction< ModelVector< Natural > > {

    public:
        Func_foldSFS(void);

        // Basic utility functions
        Func_foldSFS*                                                         clone(void) const;            //!< Clone the object
        static const std::string&                                             getClassType(void);           //!< Get Rev type
        static const TypeSpec&                                                getClassTypeSpec(void);       //!< Get class type spec
        std::string                                                           getFunctionName(void) const;  //!< Get the primary name of the function in Rev
        const TypeSpec&                                                       getTypeSpec(void) const;      //!< Get the type spec of the instance

        // Function functions you have to override
        RevBayesCore::TypedFunction< RevBayesCore::RbVector<std::int64_t> >*  createFunction(void) const;   //!< Create internal function object
        const ArgumentRules&                                                  getArgumentRules(void) const; //!< Get argument rules
    };

}

#endif
