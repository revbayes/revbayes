#ifndef RlSimplex_H
#define RlSimplex_H


#include "Probability.h"
#include "Simplex.h"
#include "RealPos.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "RlContainer.h"

#include <iostream>
#include <vector>

namespace RevLanguage {
    
    /**
     * @brief Rev simplex variable class
     *
     * This class is derived from ModelVector<RealPos> and inherits all
     * its indexing and type conversion functionality. The only changes
     * we need to introduce here are to override some vector functions
     * that are not applicable to simplices (sort, unique) and to override
     * the two functions providing assignment or external modification
     * abilities to the parser (findOrCreateElement, getElement).
     */
    
    class Simplex : public ModelObject<RevBayesCore::Simplex>, public Container {
        
    public:
        Simplex(void);                                                                                                  //!< Construct empty simplex
        Simplex(const RevBayesCore::Simplex& v);                                                                        //!< Construct simplex from double (real) vector
        Simplex(RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* c);                                                  //!< Construct simplex from DAG node

        virtual                                    ~Simplex(void);                                                      //!< Destructor        

        // Basic utility functions you have to override
        virtual RevObject*                          convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        virtual Simplex*                            clone(void) const;                                                  //!< Clone object
        RevPtr<RevVariable>                         executeMethod( std::string const &name, const std::vector<Argument> &args, bool &found );
        static const std::string&                   getClassType(void);                                                 //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                             //!< Get class type spec
        virtual const TypeSpec&                     getTypeSpec(void) const;                                            //!< Get language type of the object
        virtual Probability*                        getElement(size_t idx) const;                                       //!< Get element variable (vector of indices)
        double                                      isConvertibleTo( const TypeSpec& type, bool once ) const;
        void                                        initMethods(void);
        void                                        printValue( std::ostream& o, bool user ) const;
        virtual void                                push_back(const double &x);                                      //!< Append element to end
        virtual void                                push_back(const Probability &x);                                      //!< Append element to end
        virtual void                                push_back(const RevObject &x);                                      //!< Append element to end
        virtual size_t                              size(void) const;                                                   //!< Get number of elements in container

        // ModelVector functions that we override here to stop inappropriate actions
        void                                        sort(void);                                                         //!< Sort vector
        void                                        unique(void);                                                       //!< Remove consecutive duplicates
        
        virtual bool                                allowsModificationToCompositeContainer(void) const { return false;} //!< Does an object of this type allow transformation into a composite container?
        
        
    };
    
}


#endif

