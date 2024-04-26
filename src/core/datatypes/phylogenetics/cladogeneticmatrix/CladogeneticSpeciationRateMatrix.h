//
//  CladogeneticSpeciationRateMatrix.h
//
//  Created by Will Freyman on 8/1/17
//

#ifndef CladogeneticSpeciationRateMatrix_h
#define CladogeneticSpeciationRateMatrix_h

#include <cstddef>
#include <map>
#include <vector>
#include <iosfwd>

#include "CladogeneticProbabilityMatrix.h"
#include "Cloneable.h"
#include "DagNode.h"
#include "MemberObject.h"
#include "Printable.h"
#include "RbVector.h"
#include "Serializable.h"

namespace RevBayesCore {
    
    class CladogeneticSpeciationRateMatrix : public Cloneable, public Printable, public Serializable, public MemberObject<double>, public MemberObject<RbVector<double> >, public MemberObject<RbVector<RbVector<long> > >, public MemberObject<CladogeneticProbabilityMatrix> {
        
    public:

        CladogeneticSpeciationRateMatrix(void);
        CladogeneticSpeciationRateMatrix(size_t n);
        virtual                                                 ~CladogeneticSpeciationRateMatrix(void);
        
        bool                                                    operator==(const CladogeneticSpeciationRateMatrix &rm) const { return this == &rm; }
        bool                                                    operator!=(const CladogeneticSpeciationRateMatrix &rm) const { return !operator==(rm); }
        bool                                                    operator<(const CladogeneticSpeciationRateMatrix &rm) const { return this < &rm; }
        bool                                                    operator<=(const CladogeneticSpeciationRateMatrix &rm) const { return operator<(rm) || operator==(rm); }
        
        // pure virtual methods
        virtual CladogeneticSpeciationRateMatrix*               clone(void) const;
        virtual void                                            initFromString( const std::string &s );
        
        void                                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, CladogeneticProbabilityMatrix &rv) const;     //!< Map the member methods to internal function calls
        void                                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<long> > &rv) const;     //!< Map the member methods to internal function calls
        void                                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                                    executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;     //!< Map the member methods to internal function calls
        
        
        // virtual methods that may need to overwritten
        virtual void                                            update(void) {};
        virtual std::map<std::vector<unsigned>, double>         getEventMap(double t=0.0);
        virtual const std::map<std::vector<unsigned>, double>&  getEventMap(double t=0.0) const;
        void                                                    setEventMap(std::map<std::vector<unsigned>, double> m);
        
        // MJL: for TensorPhylo parameterization
        virtual std::vector<double>                             getSpeciationRateSumPerState(void);
        virtual const std::vector<double>&                      getSpeciationRateSumPerState(void) const;
        void                                                    setSpeciationRateSumPerState(std::vector<double> r);
        CladogeneticProbabilityMatrix                           getCladogeneticProbabilityMatrix(void);
        const CladogeneticProbabilityMatrix&                    getCladogeneticProbabilityMatrix(void) const;
        void                                                    setCladogeneticProbabilityMatrix(CladogeneticProbabilityMatrix p);
        
        
        // public methods
        size_t                                                  getNumberOfStates(void) const;      //!< Return the number of states
        size_t                                                  size(void) const;                   //!< Get the size of the rate matrix, same as the number of states
        
	json                                                    toJSON() const;
        virtual void                                            printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const;
        virtual void                                            printForSimpleStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten ) const;
        virtual void                                            printForComplexStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten ) const;
        
    protected:
        
        // protected members available for derived classes
        size_t                                                  num_states;                         //!< The number of character states
        std::map<std::vector<unsigned>, double>                 event_map;
        
        // MJL: TensorPhylo parameterization
        std::vector<double>                                     speciation_rate_sum_per_state;
        CladogeneticProbabilityMatrix                           cladogenetic_probability_matrix;
        
    };
    
};

#endif /* CladogeneticSpeciationRateMatrix_h */
