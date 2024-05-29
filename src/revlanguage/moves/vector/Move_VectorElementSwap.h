#ifndef Move_VectorElementSwap_H
#define Move_VectorElementSwap_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of a reversible jump move between a constant value and
     *  a random distributedly value.
     *
     * This class is the RevLanguage wrapper of ReversibleJumpMixtureProposal.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2014-08-18, version 1.0
     */
    template <class rlValueType>
    class Move_VectorElementSwap : public Move {
        
    public:
        
        Move_VectorElementSwap(void);                                                                                                            //!< Default constructor
        
        // Basic utility functions
        virtual Move_VectorElementSwap*             clone(void) const;                                                                              //!< Clone the object
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
        RevPtr<const RevVariable>                   neighbor;                                                                                              //!< The variable holding the real valued vector.

    };
    
}


#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "VectorElementSwapProposal.h"
#include "RbException.h"
#include "TypeSpec.h"


using namespace RevLanguage;


template <class rlValueType>
RevLanguage::Move_VectorElementSwap<rlValueType>::Move_VectorElementSwap() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <class rlValueType>
RevLanguage::Move_VectorElementSwap<rlValueType>* Move_VectorElementSwap<rlValueType>::clone(void) const
{
    
    return new Move_VectorElementSwap<rlValueType>(*this);
}


template <class rlValueType>
void RevLanguage::Move_VectorElementSwap<rlValueType>::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // only neighbors?
    bool nei = static_cast<const RlBoolean &>( neighbor->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<typename rlValueType::valueType> >* tmp = static_cast<const ModelVector<rlValueType> &>( x->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> par = tmp->getParents();
    std::vector< RevBayesCore::StochasticNode<typename rlValueType::valueType> *> n;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = par.begin(); it != par.end(); ++it)
    {
        const RevBayesCore::StochasticNode<typename rlValueType::valueType> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<typename rlValueType::valueType>* >( *it );
        if ( the_node != NULL )
        {
            n.push_back( const_cast< RevBayesCore::StochasticNode<typename rlValueType::valueType>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvVectorElementSwap because the node isn't a vector of stochastic nodes.");
        }
    }

    // now allocate a new vector-scale move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    RevBayesCore::VectorElementSwapProposal< typename rlValueType::valueType > *p = new RevBayesCore::VectorElementSwapProposal< typename rlValueType::valueType >(n, nei);
    value = new RevBayesCore::MetropolisHastingsMove(p,w);

}


/** Get Rev type of object */
template <class rlValueType>
const std::string& RevLanguage::Move_VectorElementSwap<rlValueType>::getClassType(void)
{
    
    static std::string rev_type = "Move_VectorElementSwap__" + rlValueType::getClassType();
    
    return rev_type;
}

/** Get class type spec describing type of object */
template <class rlValueType>
const RevLanguage::TypeSpec& RevLanguage::Move_VectorElementSwap<rlValueType>::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
template <class rlValueType>
std::string RevLanguage::Move_VectorElementSwap<rlValueType>::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VectorElementSwap";
    
    return c_name;
}


/** Return member rules (no members) */
template <class rlValueType>
const RevLanguage::MemberRules& RevLanguage::Move_VectorElementSwap<rlValueType>::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        move_member_rules.push_back( new ArgumentRule( "x", ModelVector<rlValueType>::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        move_member_rules.push_back( new ArgumentRule( "neighborsOnly", RlBoolean::getClassTypeSpec(), "Should we switch only neighbors or two random elements?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
template <class rlValueType>
const RevLanguage::TypeSpec& RevLanguage::Move_VectorElementSwap<rlValueType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
template <class rlValueType>
void RevLanguage::Move_VectorElementSwap<rlValueType>::printValue(std::ostream &o) const
{
    
    o << "Move_VectorElementSwap(";
    if (x != NULL)
    {
        o << x->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
template <class rlValueType>
void RevLanguage::Move_VectorElementSwap<rlValueType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "neighborsOnly" )
    {
        neighbor = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}


#endif
