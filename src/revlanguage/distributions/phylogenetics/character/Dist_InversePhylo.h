#ifndef Dist_InversePhylo_H
#define Dist_InversePhylo_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TypedDistribution.h"

namespace RevLanguage {

    template<typename valType> // Expecting AbstractHomologousDiscreteCharacterData, // could extend to other character classes
    class Dist_InversePhylo : public TypedDistribution< valType > {

    public:
        Dist_InversePhylo( void ) : TypedDistribution< valType >(), base_distribution( NULL ) {}
        
        virtual ~Dist_InversePhylo() {}

        // Basic utility functions
        // Clone the object
        Dist_InversePhylo* clone(void) const {
            return new Dist_InversePhylo(*this);
        }

        // Get Rev type
        static const std::string& getClassType(void) {
            static std::string rev_type = "Dist_InversePhylo";
            return rev_type;
        }
        
        // Get class type spec
        static const TypeSpec& getClassTypeSpec(void) {
            static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< valType >::getClassTypeSpec() ) );
            return rev_type_spec;
        }

        // Get the alternative Rev names (aliases) for the constructor function
        std::vector<std::string> getDistributionFunctionAliases( void ) const {
            // create alternative constructor function names variable that is the same for all instance of this class
            std::vector<std::string> a_names;
            a_names.push_back( "inv" );
            
            return a_names;
        }

        // Get the Rev-name for this distribution
        std::string getDistributionFunctionName(void) const {
            return "inverse";
        }  

        // Get the type spec of the instance
        const TypeSpec& getTypeSpec(void) const {
            static TypeSpec type_spec = getClassTypeSpec();
            return type_spec;
        }  

        // Get member rules (const)
        const MemberRules& getParameterRules(void) const {
            static MemberRules dist_member_rules;
            static bool rules_set = false;
            
            if ( !rules_set ) {
                dist_member_rules.push_back( new ArgumentRule( "distribution", 
                    TypedDistribution<valType>::getClassTypeSpec(), 
                    "The distribution to invert.", 
                    ArgumentRule::BY_CONSTANT_REFERENCE, 
                    ArgumentRule::ANY ) );
                rules_set = true;
            }
            
            return dist_member_rules;
        }

        // Distribution functions
        RevBayesCore::TypedDistribution<typename valType::valueType>* createDistribution(void) const;

        // Basic utility functions
        // MethodTable                                             getDistributionMethods(void) const;
        // void                                                    printValue(std::ostream& o) const;
        // void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

    protected:
        double                                                  calcLnProbability(void) const;
        void setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
            if ( name == "distribution" ) {
                base_distribution = var;
            }
            else {
                TypedDistribution< valType >::setConstParameter(name, var);
            }
        }

    private:
        RevPtr<const RevVariable> base_distribution;
    };

}

#endif