#include "AbstractCharacterData.h"

#include <string>
#include <cstddef>
#include <iostream>
#include <utility>
#include <sstream> // IWYU pragma: keep

#include "AbstractTaxonData.h"
#include "RbException.h"

using namespace RevBayesCore;

/**
 * Default constructor,
 * Does nothing except instanciating the object.
 */
AbstractCharacterData::AbstractCharacterData() {
    
}


/**
 * Copy constructor,
 * Copies the members.
 */
AbstractCharacterData::AbstractCharacterData(const AbstractCharacterData &d) :
    deletedTaxa(d.deletedTaxa),
    filename(d.filename),
    taxa(d.taxa),
    taxonMap(),
    homeologMap() {
    
    for (std::map<std::string, AbstractTaxonData*>::const_iterator it = d.taxonMap.begin(); it != d.taxonMap.end(); ++it)
    {
        const std::string &name = it->first;
        
        // add the sequence also as a member so that we can access it by name
        taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( name, it->second->clone() ) );
    }
    
    for (std::map<std::string, std::string>::const_iterator it = d.homeologMap.begin(); it != d.homeologMap.end(); ++it)
    {
        const std::string &name1 = it->first;
        const std::string &name2 = it->second;
        
        // add the sequence also as a member so that we can access it by name
        homeologMap.insert( std::pair<std::string, std::string >( name1, name2 ) );
    }
}

/**
 * Destructor.
 * Clears the data.
 */
AbstractCharacterData::~AbstractCharacterData()
{

    clear();
}



/**
 * Index (const) operator to access a TaxonData object at position i.
 *
 * \param[in]    i    The position of the TaxonData object.
 *
 * \return            The TaxonData object at position i.
 */
const AbstractTaxonData& AbstractCharacterData::operator[]( const size_t i ) const
{
    
    return getTaxonData( i );
}


AbstractCharacterData& AbstractCharacterData::operator=( const AbstractCharacterData &d )
{
    
    if ( &d != this )
    {
        clear();
        
        deletedTaxa = d.deletedTaxa;
        filename    = d.filename;
        taxa        = d.taxa;
                 
        for (std::map<std::string, AbstractTaxonData*>::const_iterator it = d.taxonMap.begin(); it != d.taxonMap.end(); ++it)
        {
            AbstractTaxonData *new_taxa_data = it->second;
            const std::string &name = it->first;
            
            // add the sequence also as a member so that we can access it by name
            taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( name, new_taxa_data->clone() ) );
        }
    }
    
    return *this;
}


///**
// * Add another character data object to this character data object.
// *
// * \param[in]    obsd    The CharacterData object that should be added.
// */
//AbstractCharacterData& AbstractCharacterData::concatenate(const AbstractCharacterData &obsd)
//{
//    
//    // check if both have the same number of taxa
//    if ( sequenceNames.size() != obsd.getNumberOfTaxa() )
//    {
//        throw RbException("Cannot add two character data objects with different number of taxa!");
//    }
//    
//    std::vector<bool> used = std::vector<bool>(obsd.getNumberOfTaxa(),false);
//    for (size_t i=0; i<sequenceNames.size(); i++ )
//    {
//        const std::string &n = sequenceNames[i];
//        AbstractTaxonData& taxon = getTaxonData( n );
//        
//        size_t idx = obsd.getIndexOfTaxon( n );
//        if ( idx != RbConstants::Size_t::inf)
//        {
//            used[idx] = true;
//            taxon.concatenate( obsd.getTaxonData( n ) );
//            
//        }
//        else
//        {
//            throw RbException() << "Cannot add two character data objects because second character data object has no taxon with name '" << n << "n'!";
//        }
//    }
//    
//    for (size_t i=0; i<used.size(); i++)
//    {
//        if ( used[i] == false )
//        {
//            throw RbException() << "Cannot add two character data objects because first character data object has no taxon with name '" << obsd.getTaxonNameWithIndex(i) << "n'!";
//        }
//    }
//    
//    return *this;
//}


/**
 * Add a sequence (TaxonData) to the character data object.
 *
 * \param[in]    n    The TaxonData object that should be added.
 */
void AbstractCharacterData::addMissingTaxon(const std::string &n) {
    
    if ( indexOfTaxonWithName(n) == -1 )
    {
        // add the sequence name to the list
        taxa.push_back( Taxon(n) );
    
        AbstractTaxonData *taxon_data = taxonMap.begin()->second->clone();
        taxon_data->setAllCharactersMissing();
    
        // add the sequence also as a member so that we can access it by name
        taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( n, taxon_data ) );
    }
}


/**
 * Add a sequence (TaxonData) to the character data object.
 *
 * \param[in]    obs    The TaxonData object that should be added.
 */
void AbstractCharacterData::addTaxonData(const AbstractTaxonData &obs)
{
	// check if the taxon is already in the data
	std::map<std::string, AbstractTaxonData* >::iterator it = taxonMap.find( obs.getTaxonName() );
	if ( it == taxonMap.end() )
	{
	    // add the sequence name to the list
	    taxa.push_back( obs.getTaxon() );
	}
	else
	{
		// remove the previously provided taxon data, so we can replace it with new taxon data
		delete it->second;
		taxonMap.erase(it);
	}


    // add the sequence as a member so that we can access it by name
    taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( obs.getTaxonName(), obs.clone() ) );
}



/**
 * Clear the object, that is, remove all TaxonData elements.
 */
void AbstractCharacterData::clear(void) {
    
    for (std::map<std::string, AbstractTaxonData*>::iterator it = taxonMap.begin(); it != taxonMap.end(); ++it)
    {
        AbstractTaxonData* d = it->second;
        delete d;
    }
    taxa.clear();
    taxonMap.clear();
}


/**
 * Exclude a taxon.
 * We don't actually delete the taxon but instead mark it for exclusion.
 *
 * \param[in]    i    The index of the taxon that will be excluded.
 */
void AbstractCharacterData::excludeTaxon(size_t i)
{
    
    if (i >= taxonMap.size())
    {
        std::stringstream o;
        o << "Only " << taxonMap.size() << " taxa in matrix";
        throw RbException( o.str() );
    }
    deletedTaxa.insert( i );
}


/**
 * Exclude a taxon.
 * We don't actually delete the taxon but instead mark it for exclusion.
 *
 * \param[in]    s    The name of the taxon that will be excluded.
 */
void AbstractCharacterData::excludeTaxon(const std::string& s) {
    
    for (size_t i = 0; i < getNumberOfTaxa(); i++)
    {
        if (s == taxa[i].getName() )
        {
            deletedTaxa.insert( i );
            break;
        }
    }
}



/**
 * Delete a taxon.
 * Remove taxon object and free up its memory.
 *
 * \param[in]    i    The index of the taxon that will be excluded.
 */
void AbstractCharacterData::deleteTaxon(size_t i) {
    
    if (i >= taxonMap.size())
    {
        std::stringstream o;
        o << "Only " << taxonMap.size() << " taxa in matrix";
        throw RbException( o.str() );
    }
    std::string name = taxa[i].getName();
    delete taxonMap[name];
    taxonMap.erase(name);
    taxa.erase(taxa.begin() + i);
}


/**
 * Delete a taxon.
 * Remove taxon object and free up its memory.
 *
 * \param[in]    s    The index of the taxon that will be excluded.
 */
void AbstractCharacterData::deleteTaxon(const std::string& s) {
    
    delete taxonMap[s];
    taxonMap.erase(s);

    std::vector<Taxon>::iterator it;
    for (it = taxa.begin(); it != taxa.end(); it++)
    {
        if (it->getName() == s) {
            break;
        }
    }
    
    taxa.erase(it);
}


/**
 * Get the file name from which the character data object was read in.
 *
 * \return    The original file name.
 */
const path& AbstractCharacterData::getFilename(void) const
{
    return filename;
}

/**
 * Get the included taxa currently stored in this object.
 *
 * \return       The included taxa.
 */
std::vector<Taxon> AbstractCharacterData::getIncludedTaxa(void) const
{
    std::vector<Taxon> included_taxa;
    size_t n_taxa_total = getNumberOfTaxa();
    
    for (size_t i=0; i<n_taxa_total; ++i)
    {
        std::set<size_t>::const_iterator pos = deletedTaxa.find( i );
        if ( pos == deletedTaxa.end() )
        {
            included_taxa.push_back( taxa[i] );
        }
    }
    
    return included_taxa;
}


/**
 * Get the index of the taxon with name 'n'.
 *
 * \par[in]     The name.
 * \return      The index.
 *
 */
size_t AbstractCharacterData::getIndexOfTaxon(const std::string &n) const
{
    
    for (size_t i=0; i<taxa.size(); ++i)
    {
        if ( taxa[i].getName() == n )
        {
            return i;
        }
    }

    throw RbException() << "Cannot find taxon '" << n << "' in the CharacterData matrix.";
}


/**
* Get a string with the JSON representation of the character data.
*
* \return     A JSON-formatted string
*/
const std::string AbstractCharacterData::getJsonRepresentation(void) const {

    //std::set<size_t>                            deletedTaxa;                                                                //!< Set of deleted taxa
    //std::string                                 fileName;                                                                   //!< The path/filename from where this matrix originated
    //std::string                                 filePath;                                                                   //!< The path/filename from where this matrix originated
    //std::map<std::string, std::string >         homeologMap;
    //std::vector<Taxon>                          taxa;                                                                       //!< names of the sequences
    //std::map<std::string, AbstractTaxonData* >  taxonMap;

    std::string jsonStr = "{\"CharacterDataMatrix\": {\n";
    
    jsonStr += std::string("   \"filePath\": \"")    + filename.parent_path().string() + std::string("\",") + '\n';
    jsonStr += std::string("   \"fileName\": \"")    + filename.parent_path().string()+ std::string("\",") + '\n';
    jsonStr += std::string("   \"dataType\": \"")    + getDataType() + "\",\n";
    jsonStr += std::string("   \"stateLabels\": \"") + getStateLabels() + "\",\n";
    
    jsonStr += "   \"taxa\": [";
    for (int i=0; i<taxa.size(); i++)
        {
        jsonStr += "\"" + taxa[i].getName() + "\"";
        if (i+1 < taxa.size())
            jsonStr += ",";
        }
    jsonStr += "],\n";

    jsonStr += "   \"deletedTaxa\": [";
    for (std::set<size_t>::iterator it = deletedTaxa.begin(); it != deletedTaxa.end(); it++)
        {
        jsonStr += *it;
        if (it != deletedTaxa.end())
            jsonStr += ",";
        }
    jsonStr += "],\n";
    
    jsonStr += "   \"data\": [\n";
    for (std::map<std::string,AbstractTaxonData* >::const_iterator it = taxonMap.begin(); it != taxonMap.end(); it++)
        {
        if (it != taxonMap.begin())
            jsonStr += ",\n";
        jsonStr += "   " + it->second->getJsonRepresentation();
        }
    jsonStr += "   ]";
    jsonStr += "}}";

    return jsonStr;
}


/**
 * Get the number of taxa currently stored in this object.
 *
 * \return       The number of taxa.
 */
size_t AbstractCharacterData::getNumberOfTaxa(void) const
{
    
    return taxa.size();
}

/**
 * Get the number of included taxa currently stored in this object.
 *
 * \return       The number of included taxa.
 */
size_t AbstractCharacterData::getNumberOfIncludedTaxa(void) const
{

    if (getNumberOfTaxa() > 0)
    {
        return getNumberOfTaxa() - deletedTaxa.size();
    }
    
    return 0;
}


/**
 * Get the percentage of missing characters n the sequences.
 *
 * \return    The percentage of missing characters.
 */
double AbstractCharacterData::getPercentageMissing( const std::string &n ) const {
    
    const AbstractTaxonData &td = getTaxonData(n);
    return td.getPercentageMissing();
}


std::string AbstractCharacterData::getStateLabels(void) {

    if (taxonMap.size() == 0)
    {
        return "";
    }
    const std::map<std::string, AbstractTaxonData* >::const_iterator& i = taxonMap.begin();
    return i->second->getStateLabels();
}

std::string AbstractCharacterData::getStateLabels(void) const {

    if (taxonMap.size() == 0)
    {
        return "";
    }
    const std::map<std::string, AbstractTaxonData* >::const_iterator& i = taxonMap.begin();
    return i->second->getStateLabels();
}


/**
 * Get the taxon with index idx.
 *
 * \param[in]    idx    The position of the taxon.
 *
 * \return              The taxon.
 */
const Taxon& AbstractCharacterData::getTaxon(size_t idx) const {
    
    return taxa[idx];
}


/**
 * Get the taxon data object with index tn.
 *
 * \return     A const reference to the taxon data object at position tn.
 */
const AbstractTaxonData& AbstractCharacterData::getTaxonData(size_t tn) const
{
    
    if ( tn >= getNumberOfTaxa() )
    {
        throw RbException( "Taxon index out of range" );
    }
    
    const std::string& name = taxa[tn].getName();
    const std::map<std::string, AbstractTaxonData* >::const_iterator& i = taxonMap.find( name );
    
    if (i != taxonMap.end() )
    {
        return *(i->second);
    }
    else
    {
        throw RbException() << "Cannot find taxon '" << name << "' in the CharacterData matrix.";
    }
}


/**
 * Get the taxon data object at position tn.
 *
 * \return     A non-const reference to the taxon data object at position tn.
 */
AbstractTaxonData& AbstractCharacterData::getTaxonData( size_t tn ) {
    
    if ( tn >= getNumberOfTaxa() )
    {
        throw RbException( "Taxon index out of range" );
    }
    
    const std::string& name = taxa[tn].getName();
    const std::map<std::string, AbstractTaxonData* >::iterator& i = taxonMap.find( name );
    
    if (i != taxonMap.end() )
    {
        return *(i->second);
    }
    else
    {
        throw RbException() << "Cannot find taxon '" << name << "' in the CharacterData matrix.";
    }
}


/**
 * Get the taxon data object with name tn.
 *
 * \return     A const reference to the taxon data object with name tn.
 */
const AbstractTaxonData& AbstractCharacterData::getTaxonData( const std::string &tn ) const {
    
    if ( tn == "" )
    {
        throw RbException("Ambiguous taxon name.");
    }
    
    const std::map<std::string, AbstractTaxonData* >::const_iterator& i = taxonMap.find(tn);
    
    if (i != taxonMap.end() )
    {
        return *(i->second);
    }
    else
    {
        throw RbException() << "Cannot find taxon '" << tn << "' in the CharacterData matrix.";
    }
}


/**
 * Get the taxon data object with name tn.
 *
 * \return     A non-const reference to the taxon data object with name tn.
 */
AbstractTaxonData& AbstractCharacterData::getTaxonData( const std::string &tn ) {
    
    if ( tn == "" )
    {
        throw RbException("Ambiguous taxon name.");
    }
    
    const std::map<std::string, AbstractTaxonData* >::iterator& i = taxonMap.find(tn);
    
    if (i != taxonMap.end() )
    {
        return *(i->second);
    }
    else
    {
        
        throw RbException() << "Cannot find taxon '" << tn << "' in the CharacterData matrix.";
    }
}


const std::map<std::string, std::string > AbstractCharacterData::getHomeologMap()
{
    return homeologMap;
}


/**
 * Get the homeolog character data currently assigned to the tip.
 *
 * \param[in] tipName        the name currently used on the tree.
 */
const std::string AbstractCharacterData::getHomeologPhase(const std::string& tipName)
{
    return homeologMap[tipName];
}


/**
 * Get the names of all taxa.
 *
 * \return     A vector of all taxon names.
 */
const std::vector<Taxon>& AbstractCharacterData::getTaxa( void ) const
{
    
    return taxa;
}



/**
 * Get the taxon name with index idx.
 *
 * \param[in]    idx    The position of the taxon.
 *
 * \return              The name of the taxon.
 */
const std::string& AbstractCharacterData::getTaxonNameWithIndex( size_t idx ) const {
    
    return taxa[idx].getName();
}


/**
 * Include a taxon.
 * Since we didn't actually deleted the taxon but marked it for exclusion
 * we can now simply remove the flag.
 *
 * \param[in]    n    The name of the taxon that will be included.
 */
void AbstractCharacterData::includeTaxon(const std::string &n) {
    
    for (size_t i = 0; i < getNumberOfTaxa(); i++)
    {
        if (n == taxa[i].getName() )
        {
            deletedTaxa.erase( i );
            break;
        }
    }
}



/**
 * Get the index of the taxon with name s.
 *
 * \param[in]    s    The name of the taxon.
 *
 * \return            The index of the taxon.
 */
size_t AbstractCharacterData::indexOfTaxonWithName( const std::string& s ) const {
    
    // search through all names
    for (size_t i=0; i<taxa.size(); i++)
    {
        if (s == taxa[i].getName() )
        {
            return i;
        }
    }
    return -1;
}


/**
 * Is the entire sequence missing, i.e., are all character '?'?
 *
 * \return     True/False whether the sequence is missing.
 */
bool AbstractCharacterData::isSequenceMissing( const std::string &n ) const {
    
    const AbstractTaxonData &td = getTaxonData(n);
    return td.isSequenceMissing();
}


/**
 * Is the taxon excluded.
 *
 * \param[in]    i    The position of the taxon in question.
 */
bool AbstractCharacterData::isTaxonExcluded(size_t i) const {
    
    std::set<size_t>::const_iterator it = deletedTaxa.find( i );
    if ( it != deletedTaxa.end() )
        return true;
    return false;
}


/**
 * Is the taxon excluded?
 *
 * \param[in]    s    The name of the taxon in question.
 */
bool AbstractCharacterData::isTaxonExcluded(const std::string& s) const {
    
    size_t i = indexOfTaxonWithName(s);
    std::set<size_t>::const_iterator it = deletedTaxa.find( i );
    if ( it != deletedTaxa.end() )
    {
        return true;
    }
    return false;
}


/**
 * Restore a taxon. We simply do not mark the taxon as excluded anymore
 *
 *
 * \param[in]    i    The position of the taxon in question.
 */
void AbstractCharacterData::restoreTaxon(size_t i) {
    
    if ( i >= getNumberOfTaxa() )
    {
        return;
    }
    deletedTaxa.erase( i );
}


/**
 * Restore a taxon. We simply do not mark the taxon as excluded anymore.
 *
 * \param[in]    s    The name of the taxon in question.
 */
void AbstractCharacterData::restoreTaxon(const std::string& s) {
    
    size_t i = indexOfTaxonWithName( s );
    deletedTaxa.erase( i );
}


/**
 * Set the original file name for this character data object.
 *
 * \param[in]    fn    The new file name.
 */
void AbstractCharacterData::setFilename(const path& fn)
{
    filename = fn;
}


/**
 * Assign character homeolog data to a tip.
 *
 * \param[in] dataName       the name currently used in the character alignment.
 * \param[in] tipName        the name to be used on the tree.
 */
void AbstractCharacterData::setHomeologPhase(const std::string& dataName, const std::string& tipName)
{
    homeologMap[tipName] = dataName;
    AbstractTaxonData& t = getTaxonData( dataName );
    t.setTaxon( Taxon(tipName) );
    size_t numTax = taxa.size();
    for (size_t i = 0; i < numTax ; ++i)
    {
        if ( taxa[i].getName() == dataName)
        {
            taxa[i] = Taxon(tipName);
            break;
        }
    }
    taxonMap.erase( dataName );
    taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( tipName, t.clone() ) );
}


/**
 * Change the name of a taxon
 *
 * \param[in] currentName    self explanatory.
 * \param[in] newName        self explanatory.
 */
void AbstractCharacterData::setTaxonName(const std::string& currentName, const std::string& newName)
{

    AbstractTaxonData& t = getTaxonData( currentName );
    t.setTaxon( Taxon(newName) );
    size_t numTax = taxa.size();
    for (size_t i = 0; i < numTax ; ++i)
        {
        if ( taxa[i].getName() == currentName)
            {
            taxa[i] = Taxon(newName);
            break;
            }
        }
    taxonMap.erase( currentName );
    taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( newName, t.clone() ) );
}


/**
 * Change the name of a taxon
 *
 * \param[in] current_name   self explanatory.
 * \param[in] new_taxon      self explanatory.
 */
void AbstractCharacterData::setTaxonObject(const std::string& current_name, const Taxon& new_taxon)
{
    
    AbstractTaxonData& t = getTaxonData( current_name );
    t.setTaxon( new_taxon );
    size_t num_taxa = taxa.size();
    for (size_t i = 0; i < num_taxa ; ++i)
    {
        if ( taxa[i].getName() == current_name)
        {
            taxa[i] = new_taxon;
            break;
        }
    }
    taxonMap.erase( current_name );
    taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( new_taxon.getName(), t.clone() ) );
}



/**
 * Print the content of the data matrix.
 */
void AbstractCharacterData::show(std::ostream &out) const {
    
    size_t nt = this->getNumberOfTaxa();
    for (size_t i=0; i<nt; i++)
        {
        const AbstractTaxonData& taxonData = this->getTaxonData(i);
        std::string taxonName = this->getTaxonNameWithIndex(i);
        size_t nc = taxonData.getNumberOfCharacters();
        std::cout << "   " << taxonName << std::endl;
        std::cout << "   ";
        for (size_t j=0; j<nc; j++)
            {
            std::string s = taxonData.getStringRepresentation(j);
            std::cout << s << " ";
            if ( (j+1) % 100 == 0 && (j+1) != nc )
                {
                std::cout << std::endl << "   ";
                }
            }
        std::cout << std::endl;
        }
}


/**
 * Switch the currently assigned homeolog phase.
 *
 * \param[in] tipName1        self explanatory. 
 * \param[in] tipName2        self explanatory.
 */
void AbstractCharacterData::switchHomeologPhase(const std::string& tipName1, const std::string& tipName2)
{
    std::string data1 = homeologMap[tipName1];
    std::string data2 = homeologMap[tipName2];
    homeologMap[tipName1] = data2;
    homeologMap[tipName2] = data1;
    AbstractTaxonData& t1 = getTaxonData( tipName1 );
    AbstractTaxonData& t2 = getTaxonData( tipName2 );
    t1.setTaxon( Taxon(tipName2) );
    t2.setTaxon( Taxon(tipName1) );
    taxonMap.erase( tipName1 );
    taxonMap.erase( tipName2 );
    taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( tipName1, &t2 ) );
    taxonMap.insert( std::pair<std::string, AbstractTaxonData* >( tipName2, &t1 ) );
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const AbstractCharacterData& x) {
    
    std::stringstream s;
    
    // Generate nice header
    o << std::endl;
    s << x.getDataType() << " character matrix with " << x.getNumberOfTaxa() << " taxa" << std::endl;
    o << s.str();
    
    for ( size_t i = 0; i < s.str().length() - 1; ++i )
        {
        o << "=";
        }
    o << std::endl;

    o << "Origination:                   " << x.getFilename().filename() << std::endl;
    o << "Number of taxa:                " << x.getNumberOfTaxa() << std::endl;
    o << "Number of included taxa:       " << x.getNumberOfIncludedTaxa() << std::endl;
    o << "Datatype:                      " << x.getDataType() << std::endl;
    o << std::endl;
    
    return o;
}





