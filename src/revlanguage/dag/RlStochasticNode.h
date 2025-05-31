#ifndef RlStochasticNode_H
#define RlStochasticNode_H

#include "StochasticNode.h"
#include "RevMemberObject.h"
#include "RlDistribution.h"

namespace RevLanguage {
    
    template<class valueType>
    class StochasticNode : public RevBayesCore::StochasticNode<valueType>, public RevMemberObject {
        
    public:
        StochasticNode(const std::string& n, RevBayesCore::TypedDistribution<valueType>* dist, Distribution* rlDist);
        StochasticNode(const StochasticNode<valueType> &n);
        virtual                            ~StochasticNode(void);
        
        StochasticNode&                     operator=( const StochasticNode &d );
        
        
        // public methods
        StochasticNode<valueType>*          clone(void) const;                                                                              //!< Clone the node
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found);         //!< Execute member method (if applicable)
        const MethodTable&                  getMethods( void ) const;                                                                       //!< Get the member methods
        Distribution&                       getRlDistribution(void);                                                                        //!< Get the Rev distribution
        const Distribution&                 getRlDistribution(void) const;                                                                  //!< Get the Rev distribution (const)
        void                                printStructureInfo(std::ostream &o, bool verbose=false) const;                                  //!< Print information on structure

    private:
        
        Distribution*                       rlDistribution;                                                                                 //!< Rev distribution
        MethodTable                         methods;
        
        
    };
    
}

#include "RlBoolean.h"
#include "RealPos.h"
#include "ModelVector.h"
#include "RlDagMemberFunction.h"

template<class valueType>
RevLanguage::StochasticNode<valueType>::StochasticNode( const std::string& n, RevBayesCore::TypedDistribution<valueType>* dist, Distribution* rlDist ) :
    RevBayesCore::StochasticNode<valueType>( n, dist ),
    rlDistribution( rlDist )
{
    
    ArgumentRules* clampArgRules = new ArgumentRules();
    clampArgRules->push_back( new ArgumentRule("x", rlDistribution->getVariableTypeSpec(), "The observed value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "clamp",  RlUtils::Void, clampArgRules) );

    ArgumentRules* integrate_out_arg_ules = new ArgumentRules();
    integrate_out_arg_ules->push_back( new ArgumentRule("x", RlBoolean::getClassTypeSpec(), "Should we integrate over this variable?", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "integrateOut", RlUtils::Void, integrate_out_arg_ules) );
    
    ArgumentRules* redrawArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "redraw", RlUtils::Void, redrawArgRules) );
    
    ArgumentRules* prob_arg_rules = new ArgumentRules();
    this->methods.addFunction( new DagMemberFunction<RealPos>( "probability", this, prob_arg_rules) );
    
    ArgumentRules* lnprob_arg_rules = new ArgumentRules();
    this->methods.addFunction( new DagMemberFunction<Real>( "lnProbability", this, lnprob_arg_rules) );
    
    ArgumentRules* ln_mixture_prob_arg_rules = new ArgumentRules();
    this->methods.addFunction( new DagMemberFunction< ModelVector<Real> >( "lnMixtureLikelihoods", this, ln_mixture_prob_arg_rules) );
    
    ArgumentRules* mixture_prob_arg_rules = new ArgumentRules();
    this->methods.addFunction( new DagMemberFunction< ModelVector<Real> >( "MixtureLikelihoods", this, mixture_prob_arg_rules) );

    ArgumentRules* setValueArgRules = new ArgumentRules();
    setValueArgRules->push_back( new ArgumentRule("x", rlDistribution->getVariableTypeSpec(), "The value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "setValue", RlUtils::Void, setValueArgRules) );
    
    ArgumentRules* unclampArgRules = new ArgumentRules();
    this->methods.addFunction( new MemberProcedure( "unclamp", RlUtils::Void, unclampArgRules) );
    
    // add the distribution member methods
    RevMemberObject* mo = dynamic_cast<RevMemberObject*>( rlDistribution );
    if ( mo != NULL)
    {
        const MethodTable &distMethods = mo->getMethods();
        methods.insertInheritedMethods( distMethods );
    }
    
}


template<class valueType>
RevLanguage::StochasticNode<valueType>::StochasticNode( const RevLanguage::StochasticNode<valueType> &n ) :
    RevBayesCore::StochasticNode<valueType>( n ),
    rlDistribution( n.rlDistribution->clone() ),
    methods( n.methods )
{
    
}



template<class valueType>
RevLanguage::StochasticNode<valueType>& RevLanguage::StochasticNode<valueType>::operator=( const RevLanguage::StochasticNode<valueType> &n )
{

    if ( this != &n )
    {
        RevBayesCore::StochasticNode<valueType>::operator=( n );

        delete rlDistribution;
        
        rlDistribution = n.rlDistribution->clone();
        methods = n.methods;
        
    }

    
}


template<class valueType>
RevLanguage::StochasticNode<valueType>::~StochasticNode( void )
{
    
    delete rlDistribution;
    
}


template<class valueType>
RevLanguage::StochasticNode<valueType>* RevLanguage::StochasticNode<valueType>::clone( void ) const
{
    
    return new StochasticNode<valueType>( *this );
}


/* Execute calls to member methods */
template <typename valueType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::StochasticNode<valueType>::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    // execute the distribution member methods
    RevMemberObject* mo = dynamic_cast<RevMemberObject*>( rlDistribution );
    if ( mo != NULL)
    {
        RevPtr<RevVariable> ret_val = mo->executeMethod(name, args, found);
        
        if ( found == true )
        {
            return ret_val;
        }
    }
    
    std::vector<RevBayesCore::DagNode*> dist_args;
    for (size_t i = 0; i < args.size(); ++i)
    {
        try
        {
            dist_args.push_back( args[i].getVariable()->getRevObject().getDagNode() );
        } catch ( ... )
        {
            // nothing to throw, just keep going
        }
    }
    
    RevPtr<RevVariable> ret_proc_val = this->distribution->executeProcedure(name, dist_args, found);
    if ( found == true )
    {
        return ret_proc_val;
    }

    
    if (name == "clamp")
    {
        
        // we found the corresponding member method
        found = true;
        
        // get the observation
//        const RevObject &tmp_obj = args[0].getVariable()->getRevObject();
//        const ModelObject<valueType> *tmp_ptr = dynamic_cast<const ModelObject<valueType> *>( &tmp_obj );
        const valueType &observation = static_cast<const ModelObject<valueType> &>( args[0].getVariable()->getRevObject() ).getValue();
        
        // clamp
        this->clamp( RevBayesCore::Cloner<valueType, IsDerivedFrom<valueType, RevBayesCore::Cloneable>::Is >::createClone( observation ) );
        
        return NULL;
    }
    else if (name == "integrateOut")
    {
        
        // we found the corresponding member method
        found = true;
        
        // get the observation
        bool tf = static_cast<const RlBoolean &>( args[0].getVariable()->getRevObject() ).getValue();
                
        // mark this node to ignore redraws
        this->setIntegratedOut( tf );
        
        return NULL;
    }
    else if (name == "redraw")
    {
        
        // we found the corresponding member method
        found = true;
        
        // manually calling redraw allows value to be set unless the node is clamped
        if (!this->isClamped())
        {
            this->setIgnoreRedraw( false );
        }
        
        // redraw the value
        this->redraw();
        
        return NULL;
    }
    else if (name == "setValue")
    {
        
        // we found the corresponding member method
        found = true;
        
        // get the observation
        const valueType &observation = static_cast<const ModelObject<valueType> &>( args[0].getVariable()->getRevObject() ).getValue();
        
        // set value
        this->setValue( RevBayesCore::Cloner<valueType, IsDerivedFrom<valueType, RevBayesCore::Cloneable>::Is >::createClone( observation ) );
        
        // mark this node to ignore redraws
        this->setIgnoreRedraw( true );
        
        return NULL;
    }
    else if (name == "unclamp")
    {
        
        // we found the corresponding member method
        found = true;
        
        // Unclamp
        this->unclamp();
        
        return NULL;
    }
    
    found = false;
    
    return NULL;
}


/**
 * Get common member methods.
 */
template<class valueType>
const RevLanguage::MethodTable& RevLanguage::StochasticNode<valueType>::getMethods( void ) const
{
    
    // return the internal value
    return methods;
}


template<class valueType>
RevLanguage::Distribution& RevLanguage::StochasticNode<valueType>::getRlDistribution( void )
{
    
    return *rlDistribution;
}


template<class valueType>
const RevLanguage::Distribution& RevLanguage::StochasticNode<valueType>::getRlDistribution( void ) const
{
    
    return *rlDistribution;
}


/** Print struct for user */
template<class valueType>
void RevLanguage::StochasticNode<valueType>::printStructureInfo( std::ostream& o, bool verbose ) const
{
    o << "_dagType      = Stochastic node (distribution)" << std::endl;
    o << "_distribution = " << rlDistribution->getRevDeclaration() << std::endl;
    o << "_clamped      = " << ( this->clamped ? "TRUE" : "FALSE" ) << std::endl;
    o << "_lnProb       = " << const_cast< StochasticNode<valueType>* >( this )->getLnProbability() << std::endl;    

    if ( verbose == true)
    {
        o << "_stored_ln_prob = ";
        if (not this->stored_ln_prob)
            o<< "EMPTY";
        else if (not *this->stored_ln_prob)
            o<< "UNCOMPUTED";
        else
            o<<**this->stored_ln_prob;

        o<< std::endl;
    }
    o << "_parents      = ";
    this->printParents( o, 16, 70, verbose );
    o << std::endl;
    
    o << "_children     = ";
    this->printChildren( o, 16, 70, verbose );
    o << std::endl;
    
    if ( verbose == true )
    {
        o << "_dagNode      = " << this->name << " <" << this << ">" << std::endl;
        o << "_refCount     = " << this->getReferenceCount() << std::endl;
    }
    
    o << std::endl;
}


#endif

