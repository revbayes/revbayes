#ifndef Dist_multispeciesCoalescentInverseGammaPrior_H
#define Dist_multispeciesCoalescentInverseGammaPrior_H

#include "MultispeciesCoalescentInverseGammaPrior.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the MultispeciesCoalescent process
     *
     * The RevLanguage wrapper of the Multispecies Coalescent with constant population process connects
     * the variables/parameters of the process and creates the internal MultispeciesCoalescent object.
     * Please read the MultispeciesCoalescent.h for more info.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_multispeciesCoalescentInverseGammaPrior :  public  TypedDistribution< ModelVector<TimeTree> > {

    public:
        Dist_multispeciesCoalescentInverseGammaPrior( void );

        // Basic utility functions
        Dist_multispeciesCoalescentInverseGammaPrior*       clone(void) const;                                                                      //!< Clone the object
        static const std::string&                           getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                         getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                  getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::MultispeciesCoalescentInverseGammaPrior*  createDistribution(void) const;

    protected:


        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:

        RevPtr<const RevVariable>                           shape;                                                                                  //!< The shape of the inverse gamma distribution
        RevPtr<const RevVariable>                           scale;                                                                                  //!< The scale of the inverse gamma distribution
        RevPtr<const RevVariable>                           species_tree;                                                                           //!< The species tree
        RevPtr<const RevVariable>                           taxa;                                                                                   //!< The taxa
        RevPtr<const RevVariable>                           num_gene_trees;                                                                         //!< The number of genes


    };

}

#endif
