#include "RlAbstractHomologousDiscreteCharacterData.h"

#include <cstddef>
#include <iostream>

#include "RlAbstractDiscreteTaxonData.h"
#include "RlDistanceMatrix.h"
#include "ArgumentRule.h"
#include "RlMatrixReal.h"
#include "MemberProcedure.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RlBoolean.h"
#include "Probability.h"
#include "RlString.h"
#include "RlSimplex.h"
#include "RlUserInterface.h"
#include "RbBitSet.h"
#include "AbstractDiscreteTaxonData.h"
#include "HomologousDiscreteCharacterData.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "DiscreteCharacterState.h"
#include "DistanceMatrix.h"
#include "MatrixReal.h"
#include "MethodTable.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Real.h"
#include "RevMemberObject.h"
#include "RevVariable.h"
#include "RlUtils.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

using namespace RevLanguage;

AbstractHomologousDiscreteCharacterData::AbstractHomologousDiscreteCharacterData(void) : ModelObject<RevBayesCore::AbstractHomologousDiscreteCharacterData>(),
HomologousCharacterData( )
{
    initMethods();
}


AbstractHomologousDiscreteCharacterData::AbstractHomologousDiscreteCharacterData( const RevBayesCore::AbstractHomologousDiscreteCharacterData &d) : ModelObject<RevBayesCore::AbstractHomologousDiscreteCharacterData>( d.clone() ),
HomologousCharacterData( )
{
    initMethods();
}


AbstractHomologousDiscreteCharacterData::AbstractHomologousDiscreteCharacterData( RevBayesCore::AbstractHomologousDiscreteCharacterData *d) : ModelObject<RevBayesCore::AbstractHomologousDiscreteCharacterData>( d ),
HomologousCharacterData( )
{
    initMethods();
}


AbstractHomologousDiscreteCharacterData::AbstractHomologousDiscreteCharacterData( RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData> *d) : ModelObject<RevBayesCore::AbstractHomologousDiscreteCharacterData>( d ),
HomologousCharacterData( )
{
    initMethods();
}



AbstractHomologousDiscreteCharacterData::~AbstractHomologousDiscreteCharacterData()
{
}



void AbstractHomologousDiscreteCharacterData::concatenate(const RevObject &d, std::string type) const
{

    const AbstractHomologousDiscreteCharacterData* tmp = dynamic_cast<const AbstractHomologousDiscreteCharacterData*>( &d );
    if ( tmp != NULL )
    {
        concatenate( *tmp, type );
    }
    else
    {
        throw RbException() << "Cannot add an object of type '" << d.getType() << "' to a character data object.";
    }

}


void AbstractHomologousDiscreteCharacterData::concatenate(const AbstractHomologousDiscreteCharacterData &d, std::string type) const
{

    // we need to make this a constant DAG node so that we can actually modify the value
    // otherwise the value might be overwritten again, e.g., if this is a deterministic node.
    //    clone_obj->makeConstantValue();

    // now concatenate
    getDagNode()->getValue().concatenate( d.getValue(), type );

}



AbstractHomologousDiscreteCharacterData* AbstractHomologousDiscreteCharacterData::clone() const
{
    return new AbstractHomologousDiscreteCharacterData( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> AbstractHomologousDiscreteCharacterData::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    if ( name == "methods" )
    {
        found = true;

        // just print the method names (including inherited methods)
        const MethodTable &m = getMethods();
        m.printValue(std::cout, true);

        return NULL;
    }

    RevPtr<RevVariable> retVal = dynamic_cast<RevMemberObject *>( dag_node )->executeMethod(name, args, found);

    if ( found == true )
    {
        return retVal;
    }

    if ( this->getDagNode() != NULL )
    {
            // set the internal value pointer
            //        setCharacterDataObject( &this->getDagNode()->getValue() );
    }

    retVal = executeCharacterDataMethod(name, args, found, &this->getValue());

    if ( found == true )
    {
        return retVal;
    }
    else if (name == "[]")
    {
        found = true;

            // get the member with give index
        const Natural& index = static_cast<const Natural&>( args[0].getVariable()->getRevObject() );

        if (this->dag_node->getValue().getNumberOfTaxa() < (size_t)(index.getValue()) )
        {
            throw RbException("Index out of bounds in []");
        }

        const RevBayesCore::AbstractDiscreteTaxonData& element = dag_node->getValue().getTaxonData(size_t(index.getValue()) - 1);

        return new RevVariable( new AbstractDiscreteTaxonData( element.clone() ) );
    }
    else if (name == "applyMissingSitesMask")
    {
        found = true;
        
        size_t                          num_taxa        = this->dag_node->getValue().getNumberOfTaxa();
        std::vector<std::vector<bool> > mask_gap        = std::vector<std::vector<bool> >(num_taxa, std::vector<bool>());
        std::vector<std::vector<bool> > mask_missing    = std::vector<std::vector<bool> >(num_taxa, std::vector<bool>());
        
        const RevBayesCore::AbstractHomologousDiscreteCharacterData& ref = static_cast<const AbstractHomologousDiscreteCharacterData&>( args[0].getVariable()->getRevObject() ).getValue();

        ref.fillMissingSitesMask(mask_gap, mask_missing);
        this->dag_node->getValue().applyMissingSitesMask(mask_gap, mask_missing);

        
        return NULL;
    }
    else if (name == "computeMultinomialProfileLikelihood")
    {
        found = true;

        double lnl = this->dag_node->getValue().computeMultinomialProfileLikelihood();

        return new RevVariable( new Real(lnl) );
    }
    else if (name == "computeSiteFrequencySpectrum")
    {
        found = true;

        bool folded = static_cast<const RlBoolean&>( args[0].getVariable()->getRevObject() ).getValue();
        const std::string& str_ambig_treat = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();

        RevBayesCore::AbstractHomologousDiscreteCharacterData::SFS_AMBIGUITY_TREATMENT ambig_treat = RevBayesCore::AbstractHomologousDiscreteCharacterData::SFS_AMBIGUITY_TREATMENT::ANCESTRAL;
        if ( str_ambig_treat == "derived" )
        {
            ambig_treat = RevBayesCore::AbstractHomologousDiscreteCharacterData::SFS_AMBIGUITY_TREATMENT::DERIVED;
        }
        else if ( str_ambig_treat == "skip" )
        {
            ambig_treat = RevBayesCore::AbstractHomologousDiscreteCharacterData::SFS_AMBIGUITY_TREATMENT::SKIP_COLUMN;
        }
        else if ( str_ambig_treat == "rescale" )
        {
            ambig_treat = RevBayesCore::AbstractHomologousDiscreteCharacterData::SFS_AMBIGUITY_TREATMENT::RESCALE;
        }

        std::vector<long> sfs = this->dag_node->getValue().computeSiteFrequencySpectrum(folded, ambig_treat);

        return new RevVariable( new ModelVector<Natural>(sfs) );
    }
    else if (name == "computeStateFrequencies")
    {
        found = true;

        RevBayesCore::MatrixReal sf = this->dag_node->getValue().computeStateFrequencies();

        return new RevVariable( new MatrixReal(sf) );
    }
    else if ( name == "expandCharacters" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        long n = static_cast<const Natural&>( argument ).getValue();

        RevBayesCore::AbstractHomologousDiscreteCharacterData *trans_data = this->dag_node->getValue().expandCharacters( n );

        return new RevVariable( new AbstractHomologousDiscreteCharacterData(trans_data) );
    }
    else if (name == "getEmpiricalBaseFrequencies")
    {
        found = true;

        std::vector<double> ebf = this->dag_node->getValue().getEmpiricalBaseFrequencies();

        return new RevVariable( new Simplex(ebf) );
    }
    else if (name == "getInvariantSiteIndices")
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        std::vector<size_t> tmp = this->dag_node->getValue().getInvariantSiteIndices( excl );
        std::vector<long> inv_vect(begin(tmp), end(tmp));

        return new RevVariable( new ModelVector<Natural>(inv_vect) );
    }
    else if (name == "getNumInvariantSites")
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        size_t n = this->dag_node->getValue().getNumberOfInvariantSites( excl );

        return new RevVariable( new Natural(n) );
    }
    else if ( name == "getStateDescriptions" )
    {
        found = true;

        std::vector<std::string> descriptions = this->dag_node->getValue().getTaxonData(0).getCharacter(0).getStateDescriptions();

        return new RevVariable( new ModelVector<RlString>(descriptions) );
    }
    else if (name == "isHomologous")
    {
        found = true;

        bool ih = this->dag_node->getValue().isHomologyEstablished();

        return new RevVariable( new RlBoolean(ih) );
    }
    else if ( name == "maxGcContent" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        double max_gc = this->dag_node->getValue().maxGcContent( excl );

        return new RevVariable( new Probability(max_gc) );
    }
    else if ( name == "maxInvariableBlockLength" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        size_t max_inv = this->dag_node->getValue().maxInvariableBlockLength( excl );

        return new RevVariable( new Natural(max_inv) );
    }
    else if ( name == "maxPairwiseDifference" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        size_t max_pd = this->dag_node->getValue().getMaxPairwiseSequenceDifference( excl );

        return new RevVariable( new Natural(max_pd) );
    }
    else if (name == "maxStates")
    {
        found = true;

        RevBayesCore::AbstractHomologousDiscreteCharacterData &v = dag_node->getValue();
        size_t nChars = v.getNumberOfCharacters();
        size_t nTaxa = v.getNumberOfTaxa();

        int max = 0;
        for (size_t i = 0; i < nChars; i++)
        {
            for (size_t j = 0; j < nTaxa; j++)
            {
                const RevBayesCore::AbstractDiscreteTaxonData& td = v.getTaxonData(j);
                if ( td.getCharacter(i).isMissingState() == false && td.getCharacter(i).isGapState() == false)
                {
                    if (td.getCharacter(i).getNumberObservedStates() > 1)
                    {
                        const RevBayesCore::RbBitSet& state = td.getCharacter(i).getState();
                        for (size_t k = 0; k < state.size(); k++)
                        {
                            if (state.test(k) && k +1 > max)
                            {
                                max = static_cast<int>(k)+1;
                            }
                        }
                    }
                    else
                    {
                        int k = int(td.getCharacter(i).getStateIndex()) + 1;
                        if (k > max)
                        {
                            max = k;
                        }
                    }
                }
            }
        }
        return new RevVariable( new Natural(max) );
    }
    else if ( name == "maxVariableBlockLength" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        size_t max_var = this->dag_node->getValue().maxVariableBlockLength( excl );

        return new RevVariable( new Natural(max_var) );
    }
    else if ( name == "meanGcContent" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        double mean_gc = this->dag_node->getValue().meanGcContent( excl );

        return new RevVariable( new Probability(mean_gc) );
    }
    else if ( name == "meanGcContentByCodonPosition" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        size_t n = size_t( static_cast<const Natural&>( argument ).getValue() );

        const RevObject& excl_argument = args[1].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( excl_argument ).getValue();

        double mean_gc = this->dag_node->getValue().meanGcContentByCodon( n, excl );

        return new RevVariable( new Probability(mean_gc) );
    }
    else if ( name == "minGcContent" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        double min_gc = this->dag_node->getValue().minGcContent( excl );

        return new RevVariable( new Probability(min_gc) );
    }
    else if ( name == "minPairwiseDifference" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        size_t min_pd = this->dag_node->getValue().getMinPairwiseSequenceDifference( excl );

        return new RevVariable( new Natural(min_pd) );
    }
    else if ( name == "getPairwiseDifference" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        RevBayesCore::DistanceMatrix pd = this->dag_node->getValue().getPairwiseSequenceDifference( excl );

        return new RevVariable( new DistanceMatrix(pd) );
    }
    else if ( name == "numInvariableBlocks" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        size_t num_blocks = this->dag_node->getValue().numInvariableSiteBlocks( excl );

        return new RevVariable( new Natural(num_blocks) );
    }
    else if ( name == "numTaxaMissingSequence" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        double percentage = static_cast<const Probability&>( argument ).getValue();

        size_t num_taxa = this->dag_node->getValue().numberTaxaMissingSequence( percentage );

        return new RevVariable( new Natural(num_taxa) );
    }
    else if (name == "removeExcludedCharacters")
    {
        found = true;
        
        this->dag_node->getValue().removeExcludedCharacters();
        
        return NULL;
    }
    else if (name == "excludeMissingSites")
    {
        found = true;
        
        this->dag_node->getValue().excludeMissingSites();
        
        return NULL;
    }
    else if (name == "replaceRandomSitesByMissingData")
    {
        found = true;
        
        double missing = static_cast<const Probability&>( args[0].getVariable()->getRevObject() ).getValue();

        this->dag_node->getValue().replaceRandomSitesByMissingData(missing);

        
        return NULL;
    }
    else if (name == "setCodonPartition")
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        RevBayesCore::AbstractHomologousDiscreteCharacterData &v = dag_node->getValue();
        size_t nChars = v.getNumberOfCharacters();

            // e.g. data.setCodonPartition(sites=v(3))
        if ( argument.isType( Natural::getClassTypeSpec() ) )
        {
            size_t n = size_t( static_cast<const Natural&>( argument ).getValue() );
            size_t i = 0; // index of included characters
            for (size_t j = 0; j < nChars; j++)
            {
                    // only set codon partition for previously included characters
                if ( !v.isCharacterExcluded(j) )
                {
                    if ( i % 3 == (n-1) )
                    {
                        v.includeCharacter(j);
                    }
                    else
                    {
                        v.excludeCharacter(j);
                    }
                    i++;
                }
            }

        }

            // e.g. data.setCodonPartition(sites=v(1,2))
        else if ( argument.isType( ModelVector<Natural>::getClassTypeSpec() ) )
        {
            const ModelVector<Natural>& x = static_cast<const ModelVector<Natural>&>( argument );
            if (x.size() == 0)
            {
                return NULL;
            }

            size_t i = 0; // index of included characters
            for (size_t j = 0; j < nChars; j++)
            {
                    // only set codon partition for previously included characters
                if ( !v.isCharacterExcluded(j) )
                {
                    bool included_codon = false;
                    for (size_t k = 0; k < x.size(); k++)
                    {
                        size_t n = x[k];
                        if ( i % 3 == (n-1) )
                        {
                            v.includeCharacter(j);
                            included_codon = true;
                            break;
                        }
                    }
                    if ( !included_codon )
                    {
                        v.excludeCharacter(j);
                    }
                    i++;
                }
            }
        }
        return NULL;
    }
    else if (name == "setNumStatesPartition")
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        RevBayesCore::AbstractHomologousDiscreteCharacterData &v = dag_node->getValue();
        size_t nChars = v.getNumberOfCharacters();
        size_t nTaxa = v.getNumberOfTaxa();

        bool warn = false;

        // e.g. data.setNumStatesPartition(2)
        size_t n = size_t( static_cast<const Natural&>( argument ).getValue() );
        for (size_t i = 0; i < nChars; i++)
        {
            // only set state number partition for previously included characters
            if ( !v.isCharacterExcluded(i) )
            {
                RevBayesCore::RbBitSet observed(v.getNumberOfStates());
                size_t max = 0;
                for (size_t j = 0; j < nTaxa; j++)
                {
                    const RevBayesCore::AbstractDiscreteTaxonData& td = v.getTaxonData(j);
                    if ( td.getCharacter(i).isMissingState() == false && td.getCharacter(i).isGapState() == false)
                    {
                        if (td.getCharacter(i).getNumberObservedStates() > 1)
                        {
                            const RevBayesCore::RbBitSet& state = td.getCharacter(i).getState();
                            for (size_t k = 0; k < state.size(); k++)
                            {
                                if (state.test(k) )
                                {
                                    observed.set(k);

                                    if ( k > max )
                                    {
                                        max = k;
                                    }
                                }
                            }
                        }
                        else
                        {
                            size_t k = td.getCharacter(i).getStateIndex();

                            observed.set(k);

                            if (k > max)
                            {
                                max = k;
                            }
                        }
                    }
                }

                if ( observed.count() != max + 1 )
                {
                    warn = true;
                }
                
                if (max + 1 == n)
                {
                    v.includeCharacter(i);
                }
                else
                {
                    v.excludeCharacter(i);
                }
            }
        }

        if ( warn == true )
        {
            std::stringstream ss;
            ss << "NOTE: " << name << "() partitions by the maximum observed state, not the total number of observed states.";
            RBOUT(ss.str());
        }

        return NULL;
    }
    else if (name == "getNumStatesVector")
    {
        found = true;

        RevBayesCore::AbstractHomologousDiscreteCharacterData &v = dag_node->getValue();
        size_t nChars = v.getNumberOfCharacters();
        size_t nTaxa = v.getNumberOfTaxa();

        bool warn = false;

        std::vector<RevBayesCore::AbstractHomologousDiscreteCharacterData*> matVec;

        for (size_t x = 0; x < v.getNumberOfStates(); x++)
        {
            matVec.push_back(v.clone());
        }

        for (size_t i = 0; i < nChars; i++)
        {
            // only set state number partition for previously included characters
            if ( !v.isCharacterExcluded(i) )
            {
                RevBayesCore::RbBitSet observed(v.getNumberOfStates());
                size_t max = 0;
                for (size_t j = 0; j < nTaxa; j++)
                {
                    const RevBayesCore::AbstractDiscreteTaxonData& td = v.getTaxonData(j);
                    if ( td.getCharacter(i).isMissingState() == false && td.getCharacter(i).isGapState() == false)
                    {
                        if (td.getCharacter(i).getNumberObservedStates() > 1)
                        {
                            const RevBayesCore::RbBitSet& state = td.getCharacter(i).getState();
                            for (size_t k = 0; k < state.size(); k++)
                            {
                                if (state.test(k) )
                                {
                                    observed.set(k);

                                    if ( k > max )
                                    {
                                        max = k;
                                    }
                                }
                            }
                        }
                        else
                        {
                            size_t k = td.getCharacter(i).getStateIndex();

                            observed.set(k);

                            if (k > max)
                            {
                                max = k;
                            }
                        }
                    }
                }

                if ( observed.count() != max + 1 )
                {
                    warn = true;
                }

                for (size_t x = 0; x < v.getNumberOfStates(); x++)
                {
                    if ( !v.isCharacterExcluded(i) && max == x ) // only consider previously included characters
                    {
                        matVec[x]->includeCharacter(i);
                    }
                    else
                    {
                        matVec[x]->excludeCharacter(i);
                    }
                }
            }
        }

        if ( warn == true )
        {
            std::stringstream ss;
            ss << "NOTE: " << name << "() partitions by the maximum observed state, not the total number of observed states.";
            RBOUT(ss.str());
        }

        ModelVector<AbstractHomologousDiscreteCharacterData>* partition = new ModelVector<AbstractHomologousDiscreteCharacterData>();

        for( size_t i = 0; i < matVec.size(); i++ )
        {
            partition->push_back(AbstractHomologousDiscreteCharacterData(matVec[i]));
        }

        return new RevVariable (partition);
    }
    else if (name == "getObsStatesVector")
    {
        found = true;

        RevBayesCore::AbstractHomologousDiscreteCharacterData &v = dag_node->getValue();
        size_t nChars = v.getNumberOfCharacters();
        size_t nTaxa = v.getNumberOfTaxa();

        bool warn = false;

        std::vector<RevBayesCore::AbstractHomologousDiscreteCharacterData*> matVec;

        for (size_t x = 0; x < v.getNumberOfStates(); x++)
        {
            matVec.push_back(v.clone());
        }

        for (size_t i = 0; i < nChars; i++)
        {
            // only set state number partition for previously included characters
            if ( !v.isCharacterExcluded(i) )
            {
                RevBayesCore::RbBitSet observed(v.getNumberOfStates());
                size_t max = 0;
                for (size_t j = 0; j < nTaxa; j++)
                {
                    const RevBayesCore::AbstractDiscreteTaxonData& taxon_data = v.getTaxonData(j);
                    if ( taxon_data.getCharacter(i).isMissingState() == false && taxon_data.getCharacter(i).isGapState() == false)
                    {
                        // Ignore ambiguous codings:
                        // each state must be unambiguously observed to be counted.
                        if (taxon_data.getCharacter(i).getNumberObservedStates() == 1)
                        {
                            size_t k = taxon_data.getCharacter(i).getStateIndex();
                            observed.set(k);
                        }
                    }
                }

                const size_t n_obs_ij = observed.count();
                for (size_t x = 0; x < v.getNumberOfStates(); x++)
                {
                    if ( !v.isCharacterExcluded(i) && n_obs_ij == x ) // only consider previously included characters
                    {
                        matVec[x]->includeCharacter(i);
                    }
                    else
                    {
                        matVec[x]->excludeCharacter(i);
                    }
                }
            }
        }

        ModelVector<AbstractHomologousDiscreteCharacterData>* partition = new ModelVector<AbstractHomologousDiscreteCharacterData>();

        for( size_t i = 0; i < matVec.size(); i++ )
        {
            partition->push_back(AbstractHomologousDiscreteCharacterData(matVec[i]));
        }

        return new RevVariable (partition);
    }
    else if ( name == "translateCharacters" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        const std::string &type = static_cast<const RlString&>( argument ).getValue();

        RevBayesCore::AbstractHomologousDiscreteCharacterData *trans_data = this->dag_node->getValue().translateCharacters( type );

        return new RevVariable( new AbstractHomologousDiscreteCharacterData(trans_data) );
    }
    else if ( name == "varGcContent" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( argument ).getValue();

        double var_gc = this->dag_node->getValue().varGcContent( excl );

        return new RevVariable( new Probability(var_gc) );
    }
    else if ( name == "varGcContentByCodonPosition" )
    {
        found = true;

        const RevObject& argument = args[0].getVariable()->getRevObject();
        size_t n = size_t( static_cast<const Natural&>( argument ).getValue() );

        const RevObject& excl_argument = args[1].getVariable()->getRevObject();
        bool excl = static_cast<const RlBoolean&>( excl_argument ).getValue();

        double var_gc = this->dag_node->getValue().varGcContentByCodon( n, excl );

        return new RevVariable( new Probability(var_gc) );
    }

    return ModelObject<RevBayesCore::AbstractHomologousDiscreteCharacterData>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& AbstractHomologousDiscreteCharacterData::getClassType(void)
{

    static std::string rev_type = "AbstractHomologousDiscreteCharacterData";

    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& AbstractHomologousDiscreteCharacterData::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( ModelObject<RevBayesCore::AbstractHomologousDiscreteCharacterData>::getClassTypeSpec() ) );

    return rev_type_spec;
}


/** Get type spec */
const TypeSpec& AbstractHomologousDiscreteCharacterData::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


void AbstractHomologousDiscreteCharacterData::initMethods( void )
{

        // add the DAG node member methods
        // note that this is a sage case because all DAG nodes are member objects
    if ( dag_node != NULL )
    {
        const MethodTable &dagMethods = dynamic_cast<RevMemberObject*>( dag_node )->getMethods();
        methods.insertInheritedMethods( dagMethods );
    }

        // insert the character data specific methods
    MethodTable charDataMethods = getCharacterDataMethods();
    methods.insertInheritedMethods( charDataMethods );

    ArgumentRules* chartypeArgRules                         = new ArgumentRules();
    ArgumentRules* comp_site_freq_spec_arg_rules            = new ArgumentRules();
    ArgumentRules* compStateFreqArgRules                    = new ArgumentRules();
    ArgumentRules* compMultiLikeArgRules                    = new ArgumentRules();
    ArgumentRules* empiricalBaseArgRules                    = new ArgumentRules();
    ArgumentRules* expandCharactersArgRules                 = new ArgumentRules();
    ArgumentRules* getNumStatesVectorArgRules               = new ArgumentRules();
    ArgumentRules* getObsStatesVectorArgRules               = new ArgumentRules();
    ArgumentRules* getPairwiseDifferenceArgRules            = new ArgumentRules();
    ArgumentRules* getStateDescriptionsArgRules             = new ArgumentRules();
    ArgumentRules* ishomologousArgRules                     = new ArgumentRules();
    ArgumentRules* invSitesArgRules                         = new ArgumentRules();
    ArgumentRules* invSiteIndicesArgRules                   = new ArgumentRules();
    ArgumentRules* mask_missing_arg_rules                   = new ArgumentRules();
    ArgumentRules* maxGcContentArgRules                     = new ArgumentRules();
    ArgumentRules* maxInvariableBlockLengthArgRules         = new ArgumentRules();
    ArgumentRules* max_states_arg_rules                     = new ArgumentRules();
    ArgumentRules* maxPairwiseDifferenceArgRules            = new ArgumentRules();
    ArgumentRules* maxVariableBlockLengthArgRules           = new ArgumentRules();
    ArgumentRules* meanGcContentArgRules                    = new ArgumentRules();
    ArgumentRules* meanGcContentByCodonPositionArgRules     = new ArgumentRules();
    ArgumentRules* minGcContentArgRules                     = new ArgumentRules();
    ArgumentRules* minPairwiseDifferenceArgRules            = new ArgumentRules();
    ArgumentRules* numInvariableBlocksArgRules              = new ArgumentRules();
    ArgumentRules* num_taxaMissingSequenceArgRules          = new ArgumentRules();
    ArgumentRules* remove_excluded_characters_arg_rules     = new ArgumentRules();
    ArgumentRules* replace_random_sites_arg_rules            = new ArgumentRules();
    ArgumentRules* exclude_missing_sites_arg_rules           = new ArgumentRules();
    ArgumentRules* setCodonPartitionArgRules                = new ArgumentRules();
    ArgumentRules* setCodonPartitionArgRules2               = new ArgumentRules();
    ArgumentRules* setNumStatesPartitionArgRules            = new ArgumentRules();
    ArgumentRules* squareBracketArgRules                    = new ArgumentRules();
    ArgumentRules* translateCharactersArgRules              = new ArgumentRules();
    ArgumentRules* varGcContentArgRules                     = new ArgumentRules();
    ArgumentRules* varGcContentByCodonPositionArgRules      = new ArgumentRules();




    
    setCodonPartitionArgRules->push_back(       new ArgumentRule("",        Natural::getClassTypeSpec()              , "The index of the codon position.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setCodonPartitionArgRules2->push_back(      new ArgumentRule("",        ModelVector<Natural>::getClassTypeSpec() , "The indicies of the codon positions.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    setNumStatesPartitionArgRules->push_back(   new ArgumentRule("",        Natural::getClassTypeSpec()              , "The number of states in this partition.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    squareBracketArgRules->push_back(           new ArgumentRule( "index" , Natural::getClassTypeSpec()              , "The index of the taxon.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );

    std::vector<std::string> sfs_ambig_options;
    sfs_ambig_options.push_back( "ancestral" );
    sfs_ambig_options.push_back( "derived" );
    sfs_ambig_options.push_back( "skip" );
    sfs_ambig_options.push_back( "rescale" );
//    argument_rules.push_back( new OptionRule( "type", new RlString("binary"), character_options, "The type of data to be constructed." ) );
    comp_site_freq_spec_arg_rules->push_back(           new ArgumentRule( "folded"           , RlBoolean::getClassTypeSpec()          , "Should we compute the folded SFS?",                   ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    comp_site_freq_spec_arg_rules->push_back(           new OptionRule( "ambiguityTreatment", new RlString("ancestral"), sfs_ambig_options, "How should we treat ambiguous characters as derived?" ) );
//    comp_site_freq_spec_arg_rules->push_back(           new ArgumentRule( "ambigAreDerived"  , RlBoolean::getClassTypeSpec()          , "Should we treat ambiguous characters as derived?",    ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    expandCharactersArgRules->push_back(                new ArgumentRule( "factor"           , Natural::getClassTypeSpec()            , "The factor by which the state space is expanded.",    ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    invSitesArgRules->push_back(                        new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    invSiteIndicesArgRules->push_back(                  new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    mask_missing_arg_rules->push_back(       new ArgumentRule("ref",        AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The reference dataset/alignment which we use for applying the mask of missing sites.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    maxGcContentArgRules->push_back(                    new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    maxInvariableBlockLengthArgRules->push_back(        new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    maxVariableBlockLengthArgRules->push_back(          new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    maxPairwiseDifferenceArgRules->push_back(           new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    minGcContentArgRules->push_back(                    new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    minPairwiseDifferenceArgRules->push_back(           new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    getPairwiseDifferenceArgRules->push_back(           new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    meanGcContentArgRules->push_back(                   new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    meanGcContentByCodonPositionArgRules->push_back(    new ArgumentRule( "index" , Natural::getClassTypeSpec()          , "The index of the codon position.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    meanGcContentByCodonPositionArgRules->push_back(    new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    numInvariableBlocksArgRules->push_back(             new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    num_taxaMissingSequenceArgRules->push_back(         new ArgumentRule( "x" ,     Probability::getClassTypeSpec()          , "The percentage/threshold for the missing sequence.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    replace_random_sites_arg_rules->push_back(       new ArgumentRule("fraction",        Probability::getClassTypeSpec(), "The fraction of sites to remove.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    translateCharactersArgRules->push_back(             new ArgumentRule( "type" ,     RlString::getClassTypeSpec()          , "The character type into which we want to translate.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    varGcContentArgRules->push_back(                    new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );
    varGcContentByCodonPositionArgRules->push_back(     new ArgumentRule( "index" , Natural::getClassTypeSpec()          , "The index of the codon position.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    varGcContentByCodonPositionArgRules->push_back(     new ArgumentRule( "excludeAmbiguous" , RlBoolean::getClassTypeSpec()          , "Should we exclude ambiguous and missing characters?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false )  ) );


    methods.addFunction( new MemberProcedure( "applyMissingSitesMask",                  RlUtils::Void,                      mask_missing_arg_rules              ) );
    methods.addFunction( new MemberProcedure( "chartype",                               RlString::getClassTypeSpec(),       chartypeArgRules                ) );
    methods.addFunction( new MemberProcedure( "computeSiteFrequencySpectrum",           ModelVector<Natural>::getClassTypeSpec(), comp_site_freq_spec_arg_rules     ) );
    methods.addFunction( new MemberProcedure( "computeStateFrequencies",                MatrixReal::getClassTypeSpec(),     compStateFreqArgRules           ) );
    methods.addFunction( new MemberProcedure( "computeMultinomialProfileLikelihood",    Real::getClassTypeSpec(),           compMultiLikeArgRules           ) );
    methods.addFunction( new MemberProcedure( "excludeMissingSites",                     RlUtils::Void,                      exclude_missing_sites_arg_rules  ) );
    methods.addFunction( new MemberProcedure( "expandCharacters",                       AbstractHomologousDiscreteCharacterData::getClassTypeSpec(),        expandCharactersArgRules         ) );
    methods.addFunction( new MemberProcedure( "getNumStatesVector"  ,                   ModelVector<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec(), getNumStatesVectorArgRules      ) );
    methods.addFunction( new MemberProcedure( "getObsStatesVector"  ,                   ModelVector<AbstractHomologousDiscreteCharacterData>::getClassTypeSpec(), getObsStatesVectorArgRules      ) );
    methods.addFunction( new MemberProcedure( "getEmpiricalBaseFrequencies",            Simplex::getClassTypeSpec(),        empiricalBaseArgRules           ) );
    methods.addFunction( new MemberProcedure( "getInvariantSiteIndices",                ModelVector<Natural>::getClassTypeSpec(), invSiteIndicesArgRules           ) );
    methods.addFunction( new MemberProcedure( "getNumInvariantSites",                   Natural::getClassTypeSpec(),        invSitesArgRules                ) );
    methods.addFunction( new MemberProcedure( "getPairwiseDifference",                  DistanceMatrix::getClassTypeSpec(), getPairwiseDifferenceArgRules       ) );
    methods.addFunction( new MemberProcedure( "getStateDescriptions",                   ModelVector<RlString>::getClassTypeSpec(), getStateDescriptionsArgRules ) );
    methods.addFunction( new MemberProcedure( "isHomologous",                           RlBoolean::getClassTypeSpec(),      ishomologousArgRules            ) );
    methods.addFunction( new MemberProcedure( "maxGcContent",                           Probability::getClassTypeSpec(),    maxGcContentArgRules                ) );
    methods.addFunction( new MemberProcedure( "maxInvariableBlockLength",               Natural::getClassTypeSpec(),        maxInvariableBlockLengthArgRules    ) );
    methods.addFunction( new MemberProcedure( "maxPairwiseDifference",                  Natural::getClassTypeSpec(),        maxPairwiseDifferenceArgRules       ) );
    methods.addFunction( new MemberProcedure( "maxStates",                              Natural::getClassTypeSpec(),        max_states_arg_rules                ) );
    methods.addFunction( new MemberProcedure( "maxVariableBlockLength",                 Natural::getClassTypeSpec(),        maxVariableBlockLengthArgRules      ) );
    methods.addFunction( new MemberProcedure( "minGcContent",                           Probability::getClassTypeSpec(),    minGcContentArgRules                ) );
    methods.addFunction( new MemberProcedure( "minPairwiseDifference",                  Natural::getClassTypeSpec(),        minPairwiseDifferenceArgRules       ) );
    methods.addFunction( new MemberProcedure( "meanGcContent",                          Probability::getClassTypeSpec(),    meanGcContentArgRules                ) );
    methods.addFunction( new MemberProcedure( "meanGcContentByCodonPosition",           Probability::getClassTypeSpec(),    meanGcContentByCodonPositionArgRules                ) );
    methods.addFunction( new MemberProcedure( "numInvariableBlocks",                    Natural::getClassTypeSpec(),        numInvariableBlocksArgRules     ) );
    methods.addFunction( new MemberProcedure( "numTaxaMissingSequence",                 Natural::getClassTypeSpec(),        num_taxaMissingSequenceArgRules ) );
    methods.addFunction( new MemberProcedure( "removeExcludedCharacters",               RlUtils::Void,                    remove_excluded_characters_arg_rules ) );
    methods.addFunction( new MemberProcedure( "replaceRandomSitesByMissingData",        RlUtils::Void,                      replace_random_sites_arg_rules   ) );
    methods.addFunction( new MemberProcedure( "setCodonPartition",                      RlUtils::Void,                      setCodonPartitionArgRules       ) );
    methods.addFunction( new MemberProcedure( "setCodonPartition",                      RlUtils::Void,                      setCodonPartitionArgRules2      ) );
    methods.addFunction( new MemberProcedure( "setNumStatesPartition",                  RlUtils::Void,                      setNumStatesPartitionArgRules   ) );
    methods.addFunction( new MemberProcedure( "translateCharacters",                    AbstractHomologousDiscreteCharacterData::getClassTypeSpec(),        translateCharactersArgRules         ) );
    methods.addFunction( new MemberProcedure( "varGcContent",                           Probability::getClassTypeSpec(),    varGcContentArgRules                ) );
    methods.addFunction( new MemberProcedure( "varGcContentByCodonPosition",            Probability::getClassTypeSpec(),    varGcContentByCodonPositionArgRules                ) );
    methods.addFunction( new MemberProcedure( "[]",                                     AbstractDiscreteTaxonData::getClassTypeSpec(), squareBracketArgRules) );

}
