#ifndef Move_DPPAllocateAuxGibbsMove_H
#define Move_DPPAllocateAuxGibbsMove_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of a scaling move on all elements of a real valued vector.
     *
     * This class is the RevLanguage wrapper of Move_DPPAllocateAuxGibbsMove.
     *
     * @author The RevBayes Development Core Team (Tracy Heath)
     * @copyright GPL version 3
     * @since 2013-11-17, version 1.0
     */
    template <typename valType>
    class Move_DPPAllocateAuxGibbsMove : public Move {
        
    public:
        
        Move_DPPAllocateAuxGibbsMove(void);                                                                                                         //!< Default constructor
        
        // Basic utility functions
        virtual Move_DPPAllocateAuxGibbsMove*       clone(void) const;                                                                              //!< Clone the object
        void                                        constructInternalObject(void);                                                                  //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                        //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                                  //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                        //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                              //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        RevPtr<const RevVariable>                   x;                                                                                              //!< The variable holding the real valued vector.
        RevPtr<const RevVariable>                   nAux;                                                                                           //!< The variable for the tuning parameter.
        
    };
    
}

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DPPAllocateAuxGibbsMove.h"
#include "ModelVector.h"
#include "Integer.h"
#include "Natural.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

template <class valType>
Move_DPPAllocateAuxGibbsMove<valType>::Move_DPPAllocateAuxGibbsMove() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <class valType>
Move_DPPAllocateAuxGibbsMove<valType>* Move_DPPAllocateAuxGibbsMove<valType>::clone(void) const
{
    
	return new Move_DPPAllocateAuxGibbsMove<valType>(*this);
}


template <class valType>
void Move_DPPAllocateAuxGibbsMove<valType>::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-scale move
    int na = (int)static_cast<const Integer &>( nAux->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<typename valType::valueType> >* tmp = static_cast<const ModelVector<valType> &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::RbVector<typename valType::valueType> > *sn = static_cast<RevBayesCore::StochasticNode< RevBayesCore::RbVector<typename valType::valueType> > *>( tmp );
    
    value = new RevBayesCore::DPPAllocateAuxGibbsMove<typename valType::valueType>(sn, na, w);
}


/** Get Rev type of object */
template <class valType>
const std::string& Move_DPPAllocateAuxGibbsMove<valType>::getClassType(void)
{
    
    static std::string rev_type = "Move_DPPAllocateAuxGibbsMove";
    
	return rev_type; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
template <class valType>
std::string Move_DPPAllocateAuxGibbsMove<valType>::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "DPPAllocateAuxGibbs";
    
    return c_name;
}


/** Get class type spec describing type of object */
template <class valType>
const RevLanguage::TypeSpec& Move_DPPAllocateAuxGibbsMove<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}



/** Return member rules (no members) */
template <class valType>
const MemberRules& Move_DPPAllocateAuxGibbsMove<valType>::getParameterRules(void) const
{
    
    static MemberRules dppMove;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        dppMove.push_back( new ArgumentRule( "x"     , ModelVector<valType>::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        dppMove.push_back( new ArgumentRule( "numAux", Integer::getClassTypeSpec()             , "The number of auxillary categories.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Integer(4) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        dppMove.insert( dppMove.end(), inheritedRules.begin(), inheritedRules.end() ); 
        
        rules_set = true;
    }
    
    return dppMove;
}

/** Get type spec */
template <class valType>
const RevLanguage::TypeSpec& Move_DPPAllocateAuxGibbsMove<valType>::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
template <class valType>
void Move_DPPAllocateAuxGibbsMove<valType>::printValue(std::ostream &o) const {
    
    o << "Move_DPPAllocateAuxGibbsMove(";
    if (x != NULL) {
        o << x->getName();
    }
    else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
template <class valType>
void Move_DPPAllocateAuxGibbsMove<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" ) {
        x = var;
    }
    else if ( name == "numAux" ) {
        nAux = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
}


#endif
