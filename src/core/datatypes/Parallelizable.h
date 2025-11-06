#ifndef Parallelizable_H
#define Parallelizable_H

#include <cstddef>

namespace RevBayesCore {
    
    
    /**
     * Interace for Parallelizable classes.
     *
     * The Parallelizable interface provides a mechanism for code parallelization.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-05-20, version 1.0
     */
    class Parallelizable {
        
    public:
        virtual                         ~Parallelizable(void) {}
        
        size_t                          getActivePID(void) const;                                               //!< Get the ID of the active process.
        size_t                          getNumberOfProcesses(void) const;                                       //!< Get the ID of the active process.
        void                            setActivePID(size_t i, size_t n);                                       //!< Set the active process id and the number of processes.
        
    protected:
        
        // hidden constrocturs
        Parallelizable();                                                                                       //!< Constructor
        Parallelizable(const Parallelizable &n);                                                                //!< Copy Constructor
        
        Parallelizable&                 operator=(const Parallelizable &p);                                     //!< Overloaded assignment operator
        
        
        // protected methods that derived classes can overwrite
        virtual void                    setActivePIDSpecialized(size_t a, size_t n);                            //!< Set the active PID in a specialized way for derived classes.
        
        
        // protected members available for derived classes
        size_t                          active_PID;
        int                             num_processes;
        int                             pid;
        bool                            process_active;
        
        
    };
    
}

#endif

