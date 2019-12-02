#ifndef Move_DPPTableValueUpdate_H
#define Move_DPPTableValueUpdate_H

#include "RlMove.h"
#include "SimpleProposal.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of a scaling move on all elements of a real valued vector.
     *
     * This class is the RevLanguage wrapper of Move_DPPTableValueUpdate.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2013-11-17, version 1.0
     */
    template <typename valType>
    class Move_DPPTableValueUpdate : public Move {
        
    public:
        
        Move_DPPTableValueUpdate(RevBayesCore::SimpleProposal<typename valType::valueType>* tmp);                                                                                                         //!< Default constructor
        
        // Basic utility functions
        virtual Move_DPPTableValueUpdate*           clone(void) const;                                                                              //!< Clone the object
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
        RevBayesCore::SimpleProposal<typename valType::valueType>*  template_proposal;
    };
    
}

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DPPTableValueUpdate.h"
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
Move_DPPTableValueUpdate<valType>::Move_DPPTableValueUpdate( RevBayesCore::SimpleProposal<typename valType::valueType>* tmp) : Move(),
    template_proposal( tmp )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <class valType>
Move_DPPTableValueUpdate<valType>* Move_DPPTableValueUpdate<valType>::clone(void) const
{
    
    return new Move_DPPTableValueUpdate<valType>(*this);
}


template <class valType>
void Move_DPPTableValueUpdate<valType>::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new vector-scale move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<typename valType::valueType> >* tmp = static_cast<const ModelVector<valType> &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::RbVector<typename valType::valueType> > *sn = static_cast<RevBayesCore::StochasticNode< RevBayesCore::RbVector<typename valType::valueType> > *>( tmp );
    
    RevBayesCore::SimpleProposal<typename valType::valueType> *p = template_proposal->clone();
    
    value = new RevBayesCore::DPPTableValueUpdate<typename valType::valueType>(sn, p, w);
}


/** 
 * Get Rev type of object 
 */
template <class valType>
const std::string& Move_DPPTableValueUpdate<valType>::getClassType(void)
{
    
    static std::string rev_type = "Move_DPPTableValueUpdate";
    return rev_type;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
template <class valType>
std::string Move_DPPTableValueUpdate<valType>::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "DPPValue" + template_proposal->getProposalName();
    
    return c_name;
}


/** 
 * Get class type spec describing type of object 
 */
template <class valType>
const RevLanguage::TypeSpec& Move_DPPTableValueUpdate<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** 
 * Return member rules (no members) 
 */
template <class valType>
const MemberRules& Move_DPPTableValueUpdate<valType>::getParameterRules(void) const
{
    
    static MemberRules rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        rules.push_back( new ArgumentRule( "x"     , ModelVector<valType>::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        rules.insert( rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return rules;
}

/** Get type spec */
template <class valType>
const RevLanguage::TypeSpec& Move_DPPTableValueUpdate<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
template <class valType>
void Move_DPPTableValueUpdate<valType>::printValue(std::ostream &o) const
{
    
    o << "Move_DPPTableValueUpdate(";
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
void Move_DPPTableValueUpdate<valType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}


#endif
