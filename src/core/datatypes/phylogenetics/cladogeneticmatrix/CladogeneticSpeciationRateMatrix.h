//
//  CladogeneticSpeciationRateMatrix.h
//
//  Created by Will Freyman on 8/1/17
//

#ifndef CladogeneticSpeciationRateMatrix_h
#define CladogeneticSpeciationRateMatrix_h

#include "Assignable.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "Printable.h"
#include <map>
#include <vector>

namespace RevBayesCore {
    
    class CladogeneticSpeciationRateMatrix : public Cloneable, public Assignable, public Printable, public Serializable {
        
    public:

        CladogeneticSpeciationRateMatrix(void);
        CladogeneticSpeciationRateMatrix(size_t n);
        virtual                                                 ~CladogeneticSpeciationRateMatrix(void);
        
        bool                                                    operator==(const CladogeneticSpeciationRateMatrix &rm) const { return this == &rm; }
        bool                                                    operator!=(const CladogeneticSpeciationRateMatrix &rm) const { return !operator==(rm); }
        bool                                                    operator<(const CladogeneticSpeciationRateMatrix &rm) const { return this < &rm; }
        bool                                                    operator<=(const CladogeneticSpeciationRateMatrix &rm) const { return operator<(rm) || operator==(rm); }
        
        // pure virtual methods
        virtual CladogeneticSpeciationRateMatrix&               assign(const Assignable &m);
        virtual CladogeneticSpeciationRateMatrix*               clone(void) const;
        virtual void                                            initFromString( const std::string &s );
        
        // virtual methods that may need to overwritten
        virtual void                                            update(void) {};
        virtual std::map<std::vector<unsigned>, double>         getEventMap(double t=0.0);
        virtual const std::map<std::vector<unsigned>, double>&  getEventMap(double t=0.0) const;
        void                                                    setEventMap(std::map<std::vector<unsigned>, double> m);
        
        // public methods
        size_t                                                  getNumberOfStates(void) const;      //!< Return the number of states
        size_t                                                  size(void) const;                   //!< Get the size of the rate matrix, same as the number of states
        
        virtual void                                            printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const;
        virtual void                                            printForSimpleStoring( std::ostream &o, const std::string &sep, int l, bool left ) const;
        virtual void                                            printForComplexStoring( std::ostream &o, const std::string &sep, int l, bool left ) const;
        
    protected:
        
        // protected members available for derived classes
        size_t                                                  num_states;                         //!< The number of character states
        std::map<std::vector<unsigned>, double>                 event_map;
        
    };
    
};

#endif /* CladogeneticSpeciationRateMatrix_h */
