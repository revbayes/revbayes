#include <cstddef>
#include <filesystem>
#include <sstream>
#include <string>

#include "MonteCarloSampler.h"
#include "Cloneable.h"
#include "DagNode.h"
#include "Parallelizable.h"
#include "RbFileManager.h"
#include "RlUserInterface.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif

using namespace RevBayesCore;


/**
 *
 */
MonteCarloSampler::MonteCarloSampler(void) : Parallelizable(),
    generation(0)
{
    
}


/**
 *
 */
MonteCarloSampler::MonteCarloSampler(const MonteCarloSampler &m) : Cloneable(m), Parallelizable(m),
    checkpoint_file_name( m.checkpoint_file_name ),
    generation( m.generation )
{
    
}


/**
 * Destructor. Nothing to do here.
 */
MonteCarloSampler::~MonteCarloSampler(void)
{
    
}


/**
 * Get the current generation number.
 */
size_t MonteCarloSampler::getCurrentGeneration( void ) const
{
    return generation;
}


void MonteCarloSampler::setCurrentGeneration( size_t g )
{
    generation = g;
}


/**
 * Public non-virtual interface entry point: does the basic work of writing out the base checkpoint file and the *_moves file, and
 * dispatches to fullCheckpoint() so the derived class can write out additional files (e.g., Mcmc::fullCheckpoint writes *_mcmc).
 */
void MonteCarloSampler::checkpoint( void )
{
    baseCheckpoint();
    
    // dispatch to the derived class
    fullCheckpoint();
}


/**
 * Performs only the work that is always needed: writing the variable names and values to a base checkpoint file, and the
 * information about the moves to a *_moves file. It intentionally does *not* call fullCheckpoint(), so callers (e.g.
 * PowerPosteriorAnalysis::burnin()) that do not want derived-class functionality (such as an *_mcmc file recording the
 * generation counter) can opt into the base step alone.
 */
void MonteCarloSampler::baseCheckpoint( void )
{
    if ( process_active == true )
    {
        // initialize variables
        std::string separator = "\t";
        bool flatten = false;
        
        createDirectoryForFile( checkpoint_file_name );
        path tmp_checkpoint_file_name = checkpoint_file_name.parent_path() / ("." + checkpoint_file_name.filename().string() + ".tmp");
        // open the stream to the file
        std::ofstream out_stream( tmp_checkpoint_file_name.string() );
        
        // first, we find the variables -- see Mcmc::resetVariableDagNodes() for reference
        std::vector<DagNode*> variable_nodes;
        
        // we only want to have each nodes once
        // this should by default happen by here we check again
        std::set<std::string> var_names;
        const std::vector<DagNode*> &n = getModel().getDagNodes();
            
        for (auto& node: n)
        {
            if ( !node->isClamped() )
            {
                if ( node->isStochastic() && !node->isHidden() )
                {
                    const std::string &name = node->getName();
                    if ( var_names.find( name ) == var_names.end() )
                    {
                        variable_nodes.push_back( node );
                        var_names.insert( name );
                    }
                }
            }
        }
        
        // we write the names of the variables
        for (std::vector<DagNode *>::const_iterator it=variable_nodes.begin(); it!=variable_nodes.end(); ++it)
        {
            // add a separator before every new element
            if ( it != variable_nodes.begin() )
            {
                out_stream << separator;
            }
            
            const DagNode* the_node = *it;
            
            // print the header
            if (the_node->getName() != "")
            {
                the_node->printName(out_stream,separator, -1, true, flatten);
            }
            else
            {
                out_stream << "Unnamed";
            }
            
        }
        out_stream << std::endl;
        
        // second, we write the values of the variables
        for (std::vector<DagNode*>::const_iterator it = variable_nodes.begin(); it != variable_nodes.end(); ++it)
        {
            // add a separator before every new element
            if ( it != variable_nodes.begin() )
            {
                out_stream << separator;
            }
            
            // get the node
            DagNode *node = *it;
            
            // print the value
            node->printValue(out_stream, separator, -1, false, false, false, flatten);
        }
        
        // clean up
        out_stream.close();
        const bool ok = out_stream.good();
        if ( !ok )
        {
            RBOUT( "Warning: failed to write checkpoint file \"" + checkpoint_file_name.string() + "\"; keeping existing file." );
            std::error_code ec;
            std::filesystem::remove(tmp_checkpoint_file_name, ec);
        }
        else
#ifdef _WIN32
            if ( MoveFileExW(tmp_checkpoint_file_name.wstring().c_str(), checkpoint_file_name.wstring().c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) == 0 )
            {
                throw RbException() << "Could not replace checkpoint file " << checkpoint_file_name;
            }
#else
        std::filesystem::rename(tmp_checkpoint_file_name, checkpoint_file_name);
#endif
        
        
        /////////
        // Next we also write the moves information into a file
        /////////
        
        // assemble the new filename
        path moves_checkpoint_file_name = appendToStem(checkpoint_file_name, "_moves");
        path tmp_moves_checkpoint_file_name = moves_checkpoint_file_name.parent_path() / ("." + moves_checkpoint_file_name.filename().string() + ".tmp");
        // open the stream to the file
        std::ofstream out_stream_moves( tmp_moves_checkpoint_file_name.string() );
        
        // get the moves
        RbVector<Move>& moves = getMoves();
        
        for (size_t i = 0; i < moves.size(); ++i)
        {
            out_stream_moves << moves[i].getMoveName();
            out_stream_moves << "(variable="                << moves[i].getDagNodes()[0]->getName();
            out_stream_moves << ",num_tried_current="       << moves[i].getNumberTriedCurrentPeriod();
            out_stream_moves << ",num_tried_total="         << moves[i].getNumberTriedTotal();
            out_stream_moves << ",num_accepted_current="    << moves[i].getNumberAcceptedCurrentPeriod();
            out_stream_moves << ",num_accepted_total="      << moves[i].getNumberAcceptedTotal();
            out_stream_moves << ",tuning_value="            << moves[i].getMoveTuningParameter();
            out_stream_moves << ")" << std::endl;
        }
        
        // clean up
        out_stream_moves.close();
        const bool ok_moves = out_stream_moves.good();
        if ( !ok_moves )
        {
            RBOUT( "Warning: failed to write checkpoint file \"" + moves_checkpoint_file_name.string() + "\"; keeping existing file." );
            std::error_code ec;
            std::filesystem::remove(tmp_moves_checkpoint_file_name, ec);
        }
        else
#ifdef _WIN32
            if ( MoveFileExW(tmp_moves_checkpoint_file_name.wstring().c_str(), moves_checkpoint_file_name.wstring().c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) == 0 )
            {
                throw RbException() << "Could not replace checkpoint file " << moves_checkpoint_file_name;
            }
#else
        std::filesystem::rename(tmp_moves_checkpoint_file_name, moves_checkpoint_file_name);
#endif
    }
}


/**
 * Public non-virtual interface entry point: does the basic work of loading the base checkpoint file and information about
 * the moves, and dispatches to fullInitializeSamplerFromCheckpoint() so the derived class can load additional files (e.g.,
 * Mcmc::fullInitializeSamplerFromCheckpoint() reads *_mcmc and restarts file monitors).
 */
void MonteCarloSampler::initializeSamplerFromCheckpoint( void )
{
    baseInitializeSamplerFromCheckpoint();
    
    fullInitializeSamplerFromCheckpoint();
}


/**
 * Performs only the work that is always needed: parsing the base checkpoint and *_moves files. It intentionally does *not*
 * call fullInitializeSamplerFromCheckpoint(), so callers (e.g. PowerPosteriorAnalysis::burnin()) that have no monitors to
 * restart and no *_mcmc file to parse can opt into the base step alone.
 */
void MonteCarloSampler::baseInitializeSamplerFromCheckpoint( void )
{
    // Open file
    std::ifstream inFile( checkpoint_file_name.string() );
    
    if ( !inFile )
    {
        throw RbException() << "Could not open file " << checkpoint_file_name;
    }
    
    // Initialize
    std::string commandLine;
    std::string delimiter = "\t";
    std::vector<std::string> parameter_names;
    std::vector<std::string> parameter_values;
    
    // our variable to store the current line of the file
    std::string line;
    
    // Command-processing loop
    while ( inFile.good() )
    {
        // Read a line
        safeGetline( inFile, line );
        
        // skip empty lines
        if (line.length() == 0)
        {
            continue;
        }
        
        // removing comments
        if (line[0] == '#')
        {
            continue;
        }
        
        break;
    }
    
    // we assume the parameter names at the first line of the file
    StringUtilities::stringSplit(line, delimiter, parameter_names);
    
    // Read a line
    safeGetline( inFile, line );
    
    // we assume the parameter values at the second line of the file
    StringUtilities::stringSplit(line, delimiter, parameter_values);
    
    // clean up
    inFile.close();
    
    size_t n_parameters = parameter_names.size();
    std::vector<DagNode*> nodes = getModel().getDagNodes();
    
    for ( size_t i = 0; i < n_parameters; ++i )
    {
        std::string parameter_name = parameter_names[i];
        
        // iterate over all DAG nodes (variables)
        for ( size_t j = 0; j < nodes.size(); ++j )
        {
            if ( nodes[j]->getName() == parameter_name )
            {
                // set the value for the variable with the last sample in the trace
                nodes[j]->setValueFromString( parameter_values[i] );
                nodes[j]->keep();
                break;
            }
        }
    }

    // We need to touch these so that their probabilities get recomputed.
    for (auto& node: nodes)
    {
        node->touch();
    }
    
    // Next we also parse the information stored in the *_moves checkpoint file
    path moves_checkpoint_file_name = appendToStem( checkpoint_file_name, "_moves" );
    
    // Open file
    std::ifstream in_file_moves( moves_checkpoint_file_name.string() );
    
    std::string line_moves;
    std::vector<std::string> stored_move_info;
    
    // Command-processing loop
    while ( in_file_moves.good() )
    {
        // Read a line
        safeGetline( in_file_moves, line_moves );
        
        if ( line_moves != "" )
        {
            stored_move_info.push_back( line_moves );
        }
    }
    
    if ( getMoves().size() != stored_move_info.size() )
    {
        throw RbException("The number of stored moves from the checkpoint file doesn't match the number of moves for this Monte Carlo sampler.");
    }
    
    for (size_t i = 0; i < getMoves().size(); ++i)
    {
        std::vector<std::string> tokens;
        StringUtilities::stringSplit( stored_move_info[i], "(", tokens);
        
        if ( getMoves()[i].getMoveName() != tokens[0] )
        {
            throw RbException("The order of the moves from the checkpoint file does not match.");
        }
        
        std::string tmp_values = tokens[1].substr(0,tokens[1].size()-1);
        std::vector<std::string> values;
        StringUtilities::stringSplit( tmp_values, ",", values);
        
        std::vector<std::string> key_value;
        StringUtilities::stringSplit( values[0], "=", key_value);
        if ( getMoves()[i].getDagNodes()[0]->getName() != key_value[1] )
        {
            throw RbException() << "The order of the moves from the checkpoint file does not match. A move working on node '" << getMoves()[i].getDagNodes()[0]->getName() << "' received a stored counterpart working on node '" << values[0] << "'.";
        }
        
        key_value.clear();
        StringUtilities::stringSplit( values[1], "=", key_value);
        getMoves()[i].setNumberTriedCurrentPeriod( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[2], "=", key_value);
        getMoves()[i].setNumberTriedTotal( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[3], "=", key_value);
        getMoves()[i].setNumberAcceptedCurrentPeriod( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[4], "=", key_value);
        getMoves()[i].setNumberAcceptedTotal( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[5], "=", key_value);
        getMoves()[i].setMoveTuningParameter( atof(key_value[1].c_str()) );
    }

    // clean up
    in_file_moves.close();
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const MonteCarloSampler& x)
{
    o << "MonteCarloSampler";
    
    return o;
}



