#ifndef Func_getNumProcesses_H
#define Func_getNumProcesses_H

#include "Procedure.h"

namespace RevLanguage {

    /**
     * @brief Rev function to get the number of MPI processes.
     *
     * For an MPI build, this is MPI_Comm_size(MPI_COMM_WORLD).
     * For a non-MPI build, this is always 1.
     */
    class Func_getNumProcesses : public Procedure {

    public:
        Func_getNumProcesses( void );

        // Basic utility functions
        Func_getNumProcesses*                       clone(void) const;                                                          //!< Clone object
        static const std::string&                   getClassType(void);                                                         //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                     //!< Get class type spec
        std::string                                 getFunctionName(void) const;
        const TypeSpec&                             getTypeSpec(void) const;                                                    //!< Get language type of the object

        // Func_getNumProcesses functions
        const ArgumentRules&                        getArgumentRules(void) const;                                               //!< Get argument rules
        const TypeSpec&                             getReturnType(void) const;                                                  //!< Get type of return val

        RevPtr<RevVariable>                         execute(void);                                                              //!< Execute function

    protected:

    };

}

#endif
