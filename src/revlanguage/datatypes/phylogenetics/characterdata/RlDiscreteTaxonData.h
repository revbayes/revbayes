#ifndef RlDiscreteTaxonData_H
#define RlDiscreteTaxonData_H

#include "DiscreteTaxonData.h"
#include "ModelObject.h"

#include <set>
#include <string>
#include <vector>


namespace RevLanguage {
    
    
    template <class rlCharType>
    class DiscreteTaxonData : public ModelObject<RevBayesCore::DiscreteTaxonData<typename rlCharType::valueType> > {
        
    public:
        DiscreteTaxonData(void);                                                                                        //!< Constructor requires character type
        DiscreteTaxonData(RevBayesCore::DiscreteTaxonData<typename rlCharType::valueType> *v);                          //!< Constructor requires character type
        virtual ~DiscreteTaxonData(){};
        
        typedef RevBayesCore::DiscreteTaxonData<typename rlCharType::valueType> valueType;
        
        // Basic utility functions
        DiscreteTaxonData*                  clone(void) const;                                                          //!< Clone object
        static const std::string&           getClassType(void);                                                         //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                     //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                    //!< Get language type of the object
        
        // Member method inits
        virtual RevPtr<RevVariable>         executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found);  //!< Override to map member methods to internal functions
        
        std::string                         getGuiName(void) { return "Taxon Data"; }
        std::string                         getGuiUnicodeSymbol(void) { return "S"; }
        std::string                         getGuiInfo(void) { return ""; }
    };
    
}


#include "ArgumentRule.h"
#include "MemberProcedure.h"
#include "Natural.h"
#include "RlBoolean.h"


template <class rlType>
RevLanguage::DiscreteTaxonData<rlType>::DiscreteTaxonData(void) : ModelObject<RevBayesCore::DiscreteTaxonData<typename rlType::valueType> >()
{

    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", rlType::getClassTypeSpec(), squareBracketArgRules) );
    
}


template <class rlType>
RevLanguage::DiscreteTaxonData<rlType>::DiscreteTaxonData( RevBayesCore::DiscreteTaxonData<typename rlType::valueType> *v) : ModelObject<RevBayesCore::DiscreteTaxonData<typename rlType::valueType> >( v )
{
    // Add method for call "x[]" as a function
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    this->methods.addFunction( new MemberProcedure( "[]", rlType::getClassTypeSpec(), squareBracketArgRules) );
    
}



template <typename rlType>
RevLanguage::DiscreteTaxonData<rlType>* RevLanguage::DiscreteTaxonData<rlType>::clone() const
{
    return new DiscreteTaxonData<rlType>( *this );
}


/* Map calls to member methods */
template <typename rlType>
RevLanguage::RevPtr<RevLanguage::RevVariable> RevLanguage::DiscreteTaxonData<rlType>::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found) {
    
    if ( name == "[]") 
    {
        found = true;
        
        // get the member with give index
        const Natural &index = static_cast<const Natural &>( args[0].getVariable()->getRevObject() );
            
        if (this->dag_node->getValue().getNumberOfCharacters() < (size_t)(index.getValue()) )
        {
            throw RbException("Index out of bounds in []");
        }
            
        RevObject* element = new rlType( this->dag_node->getValue().getCharacter( size_t(index.getValue()) - 1) );
        return new RevVariable( element );
    } 
    
    return ModelObject<RevBayesCore::DiscreteTaxonData<typename rlType::valueType> >::executeMethod( name, args, found );
}


/* Get Rev type of object */
template <typename rlType>
const std::string& RevLanguage::DiscreteTaxonData<rlType>::getClassType(void) { 
    
    static std::string rev_type = "DiscreteTaxonData";
    
	return rev_type; 
}


/* Get class type spec describing type of object */
template <typename rlType>
const RevLanguage::TypeSpec& RevLanguage::DiscreteTaxonData<rlType>::getClassTypeSpec(void) { 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ), new TypeSpec( rlType::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
template <typename rlType>
const RevLanguage::TypeSpec& RevLanguage::DiscreteTaxonData<rlType>::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}


#endif

