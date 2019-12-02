#ifndef RlNonHomologousDiscreteCharacterData_H
#define RlNonHomologousDiscreteCharacterData_H

#include "NonHomologousDiscreteCharacterData.h"
#include "ModelObject.h"
#include "RlAbstractNonHomologousDiscreteCharacterData.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    /**
     * The RevLanguage wrapper of the NonHomologousDiscreteCharacterData class.
     *
     * The RevLanguage wrapper of the discrete character data is our
     * data object that stores discrete character data such as character matrices (alignments).
     * Please read the NonHomologousDiscreteCharacterData.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-02-16, version 1.0
     *
     */
    template <class rlCharType>
    class NonHomologousDiscreteCharacterData : public AbstractNonHomologousDiscreteCharacterData {
        
    public:
        NonHomologousDiscreteCharacterData(void);                                                                                                            //!< Default constructor
        NonHomologousDiscreteCharacterData(RevBayesCore::NonHomologousDiscreteCharacterData<typename rlCharType::valueType> *v);                                          //!< Constructor with core object
        NonHomologousDiscreteCharacterData(RevBayesCore::TypedDagNode< RevBayesCore::AbstractNonHomologousDiscreteCharacterData >*d);                                     //!< Constructor with DAG node
        
//        typedef RevBayesCore::NonHomologousDiscreteCharacterData<typename rlCharType::valueType> valueType;
        
        // Basic utility functions
        NonHomologousDiscreteCharacterData*     clone(void) const;                                                                                  //!< Clone object
        static const std::string&               getClassType(void);                                                                                 //!< Get class name
        static const TypeSpec&                  getClassTypeSpec(void);                                                                             //!< Get class type spec
        const TypeSpec&                         getTypeSpec(void) const;                                                                            //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>             executeMethod(const std::string& name, const std::vector<Argument>& args, bool &f);                 //!< Override to map member methods to internal functions
        
    private:
        
        void                                    initMethods(void);
        
    };
    
}


#include "ArgumentRule.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlDiscreteTaxonData.h"


template <class rlType>
RevLanguage::NonHomologousDiscreteCharacterData<rlType>::NonHomologousDiscreteCharacterData(void) : AbstractNonHomologousDiscreteCharacterData()
{
    
    initMethods();
    
}


template <class rlType>
RevLanguage::NonHomologousDiscreteCharacterData<rlType>::NonHomologousDiscreteCharacterData( RevBayesCore::NonHomologousDiscreteCharacterData<typename rlType::valueType> *v) :
    AbstractNonHomologousDiscreteCharacterData( v )
{
    
    initMethods();
    
}


template <class rlType>
RevLanguage::NonHomologousDiscreteCharacterData<rlType>::NonHomologousDiscreteCharacterData( RevBayesCore::TypedDagNode< RevBayesCore::AbstractNonHomologousDiscreteCharacterData > *d) :
AbstractNonHomologousDiscreteCharacterData( d )
{
    
    initMethods();
    
}



template <typename charType>
RevLanguage::NonHomologousDiscreteCharacterData<charType>* RevLanguage::NonHomologousDiscreteCharacterData<charType>::clone() const {
    
    return new NonHomologousDiscreteCharacterData<charType>( *this );
}


/* Map calls to member methods */
template <typename rlType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::NonHomologousDiscreteCharacterData<rlType>::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "[]")
    {
        found = true;
        
        // get the member with give index
        const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );
        
        if (this->dag_node->getValue().getNumberOfTaxa() < (size_t)(index.getValue()) )
        {
            throw RbException("Index out of bounds in []");
        }
        
        const RevBayesCore::DiscreteTaxonData<typename rlType::valueType>& element = static_cast< RevBayesCore::NonHomologousDiscreteCharacterData<typename rlType::valueType>& >( this->dag_node->getValue() ).getTaxonData(size_t(index.getValue()) - 1);
        
        return new RevVariable( new DiscreteTaxonData<rlType>( new RevBayesCore::DiscreteTaxonData<typename rlType::valueType>( element ) ) );
    }
    
    return AbstractNonHomologousDiscreteCharacterData::executeMethod( name, args, found );
}


/* Get class name of object */
template <typename rlType>
const std::string& RevLanguage::NonHomologousDiscreteCharacterData<rlType>::getClassType(void)
{
    
    static std::string revClassType = "NonHomologousDiscreteCharacterData<" + rlType::getClassType() + ">";
    
    return revClassType;
}

/* Get class type spec describing type of object */
template <typename rlType>
const RevLanguage::TypeSpec& RevLanguage::NonHomologousDiscreteCharacterData<rlType>::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( AbstractNonHomologousDiscreteCharacterData::getClassTypeSpec() ), new TypeSpec( rlType::getClassTypeSpec() ) );
    
    return revClassTypeSpec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
template <typename rlType>
const RevLanguage::TypeSpec& RevLanguage::NonHomologousDiscreteCharacterData<rlType>::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


/* Get class name of object */
template <typename rlType>
void RevLanguage::NonHomologousDiscreteCharacterData<rlType>::initMethods( void )
{
    
    MethodTable methods = AbstractNonHomologousDiscreteCharacterData::getCharacterDataMethods();

    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberProcedure( "[]", DiscreteTaxonData<rlType>::getClassTypeSpec(), squareBracketArgRules) );

}


#endif

