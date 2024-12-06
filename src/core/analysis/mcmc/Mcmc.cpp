#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "DagNode.h"
#include "Mcmc.h"
#include "MoveSchedule.h"
#include "RandomMoveSchedule.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RlUserInterface.h"
#include "SingleRandomMoveSchedule.h"
#include "SequentialMoveSchedule.h"
#include "AbstractFileMonitor.h"
#include "Model.h"
#include "Monitor.h"
#include "MonteCarloAnalysisOptions.h"
#include "MonteCarloSampler.h"
#include "Move.h"
#include "RbConstIterator.h"
#include "RbConstIteratorImpl.h"
#include "RbFileManager.h"
#include "RbIterator.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RbSettings.h" // for logMCMC setting
#include "StringUtilities.h"

#ifdef RB_MPI
#include <mpi.h>
#endif


using namespace RevBayesCore;


/**
 * Constructor. We create an independent copy of the model and thus of all DAG nodes.
 * Someone might have wanted to run another MCMC with different settings on the same model.
 * Thus we also create our own copies of the monitors and moves.
 *
 * \param[in]    m    The model containing all DAG nodes.
 * \param[in]    mvs  The vector of moves.
 * \param[in]    mons The vector of monitors.
 */
Mcmc::Mcmc(const Model& m, const RbVector<Move> &mvs, const RbVector<Monitor> &mons, size_t ntries) : MonteCarloSampler(),
    chain_active( true ),
    chain_likelihood_heat( 1.0 ),
    chain_posterior_heat( 1.0 ),
    chain_prior_heat( 1.0 ),
    chain_idx( 0 ),
    model( m.clone() ),
    monitors( mons ),
    moves( mvs ),
    num_init_attempts(ntries),
    schedule(NULL),
    schedule_type("random")
{
    // create an independent copy of the model, monitors and moves
    replaceDag(mvs,mons);
    
    tuningInfo ti;
    ti.num_tried_current_period = 0;
    ti.num_tried_total = 0;
    ti.num_accepted_current_period = 0;
    ti.num_accepted_total = 0;
    ti.tuning_parameter = RbConstants::Double::neginf;
    
    moves_tuningInfo = std::vector<tuningInfo> (moves.size(), ti);
    for (size_t i = 0; i < moves.size(); ++i)
    {
        moves_tuningInfo[i].num_tried_current_period = moves[i].getNumberTriedCurrentPeriod();
        moves_tuningInfo[i].num_tried_total = moves[i].getNumberTriedTotal();
        moves_tuningInfo[i].num_accepted_current_period = moves[i].getNumberAcceptedCurrentPeriod();
        moves_tuningInfo[i].num_accepted_total = moves[i].getNumberAcceptedTotal();
        moves_tuningInfo[i].tuning_parameter = moves[i].getMoveTuningParameter();
    }
    
    initializeSampler();
    initializeMonitors();

}


/**
  * Copy constructor. For more details see the constructor.
  *
  * \param[in]    m    The MonteCarloSampler object to copy.
  */
Mcmc::Mcmc(const Mcmc &m) : MonteCarloSampler(m),
    chain_active( m.chain_active ),
    chain_likelihood_heat( m.chain_likelihood_heat ),
    chain_posterior_heat( m.chain_posterior_heat ),
    chain_prior_heat( m.chain_prior_heat ),
    chain_idx( m.chain_idx ),
    model( m.model->clone() ),
    monitors( m.monitors ),
    moves( m.moves ),
    num_init_attempts( m.num_init_attempts ),
    schedule( NULL ),
    schedule_type( m.schedule_type )
{
    
    // temporary references
    const RbVector<Monitor>& mons = m.monitors;
    const RbVector<Move>& mvs = m.moves;
    
    
    // create an independent copy of the model, monitors and moves
    replaceDag(mvs,mons);
    
    moves_tuningInfo = m.moves_tuningInfo;
    
    initializeSampler();
    initializeMonitors();
}


/**
 * Destructor. Frees the DAG nodes (the model), moves, monitor and the move schedule.
 */
Mcmc::~Mcmc(void)
{
    
    
    // delete the move schedule
    delete schedule;
    
    // delete the model
    delete model;
    
}


/**
 * Copy constructor. For more details see the constructor.
 *
 * \param[in]    m    The MonteCarloSampler object to copy.
 */
Mcmc& Mcmc::operator=(const Mcmc &m)
{
    
    if ( this != &m )
    {
        delete model;
        model = m.model->clone();
        
        // temporary references
        const RbVector<Monitor>& mons = m.monitors;
        const RbVector<Move>& mvs = m.moves;
    
    
        // create an independent copy of the model, monitors and moves
        replaceDag(mvs,mons);
    
        initializeSampler();
        initializeMonitors();
    }
    
    return *this;
}


/**
 * Add an extension to the name of the monitor.
 * We tell this to all our monitors.
 */
void Mcmc::addFileMonitorExtension(const std::string &s, bool dir)
{
    
    // tell each monitor
    for (RbIterator<Monitor> it=monitors.begin(); it!=monitors.end(); ++it)
    {
        it->addFileExtension( s, dir );
    }
    
}


void Mcmc::addMonitor(const Monitor &m)
{
    
    monitors.push_back( m );
    
}


Mcmc* Mcmc::clone( void ) const
{
    
    return new Mcmc( *this );
}


void Mcmc::checkpoint( void ) const
{
    // initialize variables
    std::string separator = "\t";
    bool flatten = false;
    
    createDirectoryForFile( checkpoint_file_name );
    
    // open the stream to the file
    std::ofstream out_stream( checkpoint_file_name.string() );

    // first, we write the names of the variables
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
    
    
    /////////
    // Now we also write the MCMC information into a file
    /////////

    // assemble the new filename
    path mcmc_checkpoint_file_name = appendToStem(checkpoint_file_name, "_mcmc");
    
    // open the stream to the file
    std::ofstream out_stream_mcmc( mcmc_checkpoint_file_name.string() );
    out_stream_mcmc << "iter = " << generation << std::endl;
    
    // clean up
    out_stream_mcmc.close();
    
    
    /////////
    // Next we also write the moves information into a file
    /////////
    
    // assemble the new filename
    path moves_checkpoint_file_name = appendToStem(checkpoint_file_name, "_moves");
    
    // open the stream to the file
    std::ofstream out_stream_moves( moves_checkpoint_file_name.string() );
    
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
}


/**
 * Disable all screen monitors. This means we simply delete it.
 */
void Mcmc::disableScreenMonitor( bool all, size_t rep )
{
    
    // tell each monitor
    for (size_t i=0; i < monitors.size(); ++i)
    {
     
        if ( all == true || rep > 0 || process_active == false )
        {

            bool is = monitors[i].isScreenMonitor();
            if ( is == true )
            {
                monitors[i].disable();
            }
            
        }
        
    }
    
}


/**
 * Finish the monitors which will close the output streams.
 */
void Mcmc::finishMonitors( size_t n_reps, MonteCarloAnalysisOptions::TraceCombinationTypes tc )
{
    
    // iterate over all monitors
    for (size_t i=0; i<monitors.size(); ++i)
    {
        
        // if this chain is active, then close the stream
        if ( chain_active == true && process_active == true )
        {
            monitors[i].closeStream();
            
            // combine results if we used more than one replicate
            if ( n_reps > 1 && tc != MonteCarloAnalysisOptions::NONE )
            {
                monitors[i].combineReplicates( n_reps, tc );
            }
            
        }
        
    }
    
}


/**
 * Get the heat of the likelihood of this chain.
 */
double Mcmc::getChainLikelihoodHeat(void) const
{
    return chain_likelihood_heat;
}


/**
 * Get the heat of the posterior of this chain.
 */
double Mcmc::getChainPosteriorHeat(void) const
{
    return chain_posterior_heat;
}


/**
 * Get the heat of the prior of this chain.
 */
double Mcmc::getChainPriorHeat(void) const
{
    return chain_prior_heat;
}


/**
 * Get the index of this chain.
 */
size_t Mcmc::getChainIndex(void) const
{
    return chain_idx;
}


/**
 * Is the current chain active?
 */
bool Mcmc::isChainActive(void)
{
    return chain_active;
}


/**
 * Get the model instance.
 */
const Model& Mcmc::getModel( void ) const
{
    
    return *model;
}


/**
 * Get the joint posterior probability of the current state for this model.
 * Note that the joint posterior is the true, unscaled and unheated value.
 */
double Mcmc::getModelLnProbability(bool likelihood_only)
{
    double pp = 0.0;
    
    const std::vector<DagNode*> &n = model->getDagNodes();
    for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
    {
        
        DagNode *the_node = *it;
        if (likelihood_only == false)
        {
            pp += the_node->getLnProbability();
        }
        else if (the_node->isClamped() == true )
        {
            pp += the_node->getLnProbability();
        }

    }
    
    return pp;
}


/**
 * Get the vector of monitors for this sampler.
 */
RbVector<Monitor>& Mcmc::getMonitors(void)
{
    return monitors;
}


/**
 * Get the vector of moves for this sampler.
 */
RbVector<Move>& Mcmc::getMoves(void)
{
    return moves;
}


std::vector<Mcmc::tuningInfo> Mcmc::getMovesTuningInfo(void)
{
    
    if (moves_tuningInfo.size() != moves.size())
    {
        throw RbException( "The number of moves does not match the length of tuning information structures." );
    }
    
    // iterate over the moves
    for (size_t i = 0; i < moves.size(); ++i)
    {
        moves_tuningInfo[i].num_tried_current_period    = moves[i].getNumberTriedCurrentPeriod();
        moves_tuningInfo[i].num_tried_total             = moves[i].getNumberTriedTotal();
        moves_tuningInfo[i].num_accepted_current_period = moves[i].getNumberAcceptedCurrentPeriod();
        moves_tuningInfo[i].num_accepted_total          = moves[i].getNumberAcceptedTotal();
        
        double tmp_tuningParameter = moves[i].getMoveTuningParameter();
        
        if ((std::isnan(moves_tuningInfo[i].tuning_parameter) == true && std::isnan(tmp_tuningParameter) == false) || (std::isnan(moves_tuningInfo[i].tuning_parameter) == false && std::isnan(tmp_tuningParameter) == true))
        {
            throw RbException( "The tunability of some moves changed." );
        }
        else if (std::isnan(tmp_tuningParameter) == false)
        {
            moves_tuningInfo[i].tuning_parameter = moves[i].getMoveTuningParameter();
        }
    }
    
    return moves_tuningInfo;
}


/**
 * Get a const-reference move-schedule for this sampler.
 */
const MoveSchedule& Mcmc::getSchedule(void) const
{
    return *schedule;
}


/**
 * Get a non-const reference to the move-schedule of this sampler.
 */
MoveSchedule& Mcmc::getSchedule(void)
{
    return *schedule;
}

/**
 * Get the schedule type of this sampler.
 */
const std::string& Mcmc::getScheduleType( void ) const
{
    return schedule_type;
}


std::string Mcmc::getStrategyDescription( void ) const
{
    
    std::string description = "";
    std::stringstream stream;
    if ( schedule_type == "single" )
    {
        stream << "The simulator uses " << moves.size() << " different moves, with a single move picked randomly per iteration" << std::endl;
    }
    else if ( schedule_type == "random" )
    {
        stream << "The simulator uses " << moves.size() << " different moves in a random move schedule with " << schedule->getNumberMovesPerIteration() << " moves per iteration" << std::endl;
    }
    else if ( schedule_type == "sequential" )
    {
        stream << "The simulator uses " << moves.size() << " different moves in a sequential move schedule with " << schedule->getNumberMovesPerIteration() << " moves per iteration" << std::endl;
    }
    description = stream.str();

    return description;
}


void Mcmc::initializeSampler( bool prior_only )
{
    std::vector<DagNode *> &dag_nodes = model->getDagNodes();
    std::vector<DagNode *> ordered_stoch_nodes = model->getOrderedStochasticNodes(  );

    // Get rid of previous move schedule, if any
    if ( schedule != NULL )
    {
        delete schedule;
    }
    schedule = NULL;

    // Get initial ln_probability of model

    // first we touch all nodes so that the likelihood is dirty
    for (auto the_node: dag_nodes)
    {
        the_node->setMcmcMode( true );
        if (the_node->isClamped())
            the_node->setPriorOnly( prior_only );
        the_node->touch();
    }


    if ( chain_active == false )
    {

        for (auto the_node: ordered_stoch_nodes)
        {
            if ( the_node->isClamped() == false && the_node->isStochastic() == true )
            {

                the_node->redraw();
                the_node->reInitialized();
    
            }
            else if ( the_node->isClamped() == true )
            {
                // make sure that the clamped node also recompute their probabilities
                the_node->touch();
            }
    
        }
        
    }


    int num_tries     = 0;
    double ln_probability = 0.0;
    for ( ; num_tries < num_init_attempts; ++num_tries )
    {
        // a flag if we failed to find a valid starting value
        bool failed = false;

        ln_probability = 0.0;
        for (auto the_node: dag_nodes)
            the_node->touch();

        for (auto the_node: dag_nodes)
        {
            double ln_prob = the_node->getLnProbability();

            if ( RbMath::isAComputableNumber(ln_prob) == false )
            {
                std::stringstream ss;
                ss << "Could not compute lnProb for node '" << the_node->getName() << "': lnProb = "<< ln_prob << std::endl;
                std::ostringstream o1;
                the_node->printValue( o1, "," );
                ss << StringUtilities::oneLiner( o1.str(), 54 ) << std::endl;

                ss << std::endl;
                RBOUT( ss.str() );

                // set the flag
                failed = true;
            }
            ln_probability += ln_prob;
        }

        // now we keep all nodes so that the likelihood is stored
        for (auto the_node: dag_nodes)
        {
            the_node->keep();
        }

        if ( failed == true )
        {
            RBOUT( "Drawing new initial states ... " );
            for (auto the_node: ordered_stoch_nodes)
            {
                if ( the_node->isClamped() == false && the_node->isStochastic() == true )
                {
                    the_node->redraw();
                    the_node->reInitialized();
                    
                }
                else if ( the_node->isClamped() == true )
                {
                    // make sure that the clamped node also recompute their probabilities
                    the_node->reInitialized();
                    the_node->touch();
                }
            }

            for (auto the_node: ordered_stoch_nodes)
		if (the_node->isClamped())
		    the_node->keep();
        }
        else
        {
            break;
        }

    }
    
    if ( num_tries == num_init_attempts )
    {
        std::stringstream msg;
        msg << "Unable to find a starting state with computable probability";
        if ( num_tries > 1 )
        {
            msg << " after " << num_tries << " tries";
        }
        throw RbException( msg.str() );
        
    }
    
    // Create the move scheduler
    if ( schedule_type == "sequential" )
    {
        schedule = new SequentialMoveSchedule( &moves );
    }
    else if ( schedule_type == "single" )
    {
        schedule = new SingleRandomMoveSchedule( &moves );
    }
    else
    {
        schedule = new RandomMoveSchedule( &moves );
    }
    
    generation = 0;
    
    resetVariableDagNodes();
}


void Mcmc::initializeSamplerFromCheckpoint( void )
{
    
    //    size_t n_samples = traces[0].size();
    size_t last_generation = 0;
    //    size_t n_traces = traces.size();
    
    std::vector<std::string> parameter_names;
    std::vector<std::string> parameter_values;
    
    
    // check that the file/path name has been correctly specified
    if ( not is_regular_file( checkpoint_file_name) )
    {
        std::string errorStr = "";
        formatError( checkpoint_file_name, errorStr );
        throw RbException(errorStr);
    }
    
    // Open file
    std::ifstream inFile( checkpoint_file_name.string() );
    
    if ( !inFile )
    {
        throw RbException()<<"Could not open file "<<checkpoint_file_name;
    }
    
    // Initialize
    std::string commandLine;
    std::string delimiter = "\t";
    
    // our variable to store the current line of the file
    std::string line;
    
    // Command-processing loop
    while ( inFile.good() )
    {
        
        // Read a line
        safeGetline( inFile, line );
        
        // skip empty lines
        //line = stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
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
    for(auto& node: nodes)
    {
        node->touch();
    }

    // assemble the new filename
    path mcmc_checkpoint_file_name = appendToStem( checkpoint_file_name, "_mcmc");

    // Open file
    std::ifstream in_file_mcmc( mcmc_checkpoint_file_name.string() );

    std::string line_mcmc;
    std::map<std::string, std::string> mcmc_pars;
    // Command-processing loop
    while ( in_file_mcmc.good() )
    {
        
        // Read a line
        safeGetline( in_file_mcmc, line_mcmc );
        
        if ( line_mcmc != "" )
        {
            std::vector<std::string> key_value;
            StringUtilities::stringSplit(line_mcmc, " = ", key_value);

            mcmc_pars.insert( std::pair<std::string, std::string>(key_value[0],key_value[1]) );
        }
        
    }
    last_generation = StringUtilities::asIntegerNumber( mcmc_pars["iter"] );
    
    // clean up
    in_file_mcmc.close();
    
    
    // we also need to tell our monitors to append after the last sample
    // set iteration num
    setCurrentGeneration( last_generation );
        
    for (size_t j = 0; j < monitors.size(); ++j)
    {
        if ( monitors[j].isFileMonitor() )
        {
            // set file monitors to append
            AbstractFileMonitor* m = dynamic_cast< AbstractFileMonitor *>( &monitors[j] );
            m->setAppend(true);
        }
    }
    
    
    /////////
    // Next we also write the moves information into a file
    /////////
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
    
    if ( moves.size() != stored_move_info.size() )
    {
        throw RbException("The number of stored moves from the checkpoint file doesn't match the number of moves for this MCMC analysis.");
    }
    
    for (size_t i = 0; i < moves.size(); ++i)
    {
        std::vector<std::string> tokens;
        StringUtilities::stringSplit( stored_move_info[i], "(", tokens);
        
        if ( moves[i].getMoveName() != tokens[0] )
        {
            throw RbException("The order of the moves from the checkpoint file does not match.");
        }
        
        std::string tmp_values = tokens[1].substr(0,tokens[1].size()-1);
        std::vector<std::string> values;
        StringUtilities::stringSplit( tmp_values, ",", values);
        
        std::vector<std::string> key_value;
        StringUtilities::stringSplit( values[0], "=", key_value);
        if ( moves[i].getDagNodes()[0]->getName() != key_value[1] )
        {
            throw RbException() << "The order of the moves from the checkpoint file does not match. A move working on node '" << moves[i].getDagNodes()[0]->getName() << "' received a stored counterpart working on node '" << values[0] << "'.";
        }
        
        key_value.clear();
        StringUtilities::stringSplit( values[1], "=", key_value);
        moves[i].setNumberTriedCurrentPeriod( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[2], "=", key_value);
        moves[i].setNumberTriedTotal( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[3], "=", key_value);
        moves[i].setNumberAcceptedCurrentPeriod( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[4], "=", key_value);
        moves[i].setNumberAcceptedTotal( StringUtilities::asIntegerNumber(key_value[1]) );
        
        key_value.clear();
        StringUtilities::stringSplit( values[5], "=", key_value);
        moves[i].setMoveTuningParameter( atof(key_value[1].c_str()) );
        
    }

    // clean up
    in_file_moves.close();
    
}


void Mcmc::initializeMonitors(void)
{
    
    for (size_t i=0; i<monitors.size(); ++i)
    {
        monitors[i].setModel( model );
    }

}


void Mcmc::monitor(unsigned long g)
{
    
    if ( chain_active == true && process_active == true )
    {
        // Monitor
        for (size_t i = 0; i < monitors.size(); ++i)
        {
            
            monitors[i].monitor( g );
        
        }
        
    }
    
}


void Mcmc::nextCycle(bool advance_cycle)
{
    int logMCMC = RbSettings::userSettings().getLogMCMC();

    size_t proposals = size_t( round( schedule->getNumberMovesPerIteration() ) );
    
    for (size_t i=0; i<proposals; ++i)
    {
        
        // Get the move
        Move& the_move = schedule->nextMove( generation );

	if (logMCMC >= 1)
	{
	    std::vector<std::string> node_names;
	    for(auto node: the_move.getDagNodes())
		node_names.push_back(node->getName());

	    std::cerr<<"\ngeneration = "<<generation<<"    proposal = "<<i+1<<"/"<<proposals<<"    "<<the_move.getMoveName()<<"("<<StringUtilities::join(node_names,",")<<")\n";
	}

        // Perform the move
        the_move.performMcmcStep( chain_prior_heat, chain_likelihood_heat, chain_posterior_heat );
        
    }
    
    
    // advance gen cycle if needed (i.e. run()==true, burnin()==false)
    if ( advance_cycle == true )
    {
        ++generation;
    }

}



void Mcmc::printOperatorSummary(bool current_period)
{
    
    if ( process_active == true )
    {
        // printing the moves summary
        std::cout << std::endl;
        std::cout << "                  Name                  | Param              |  Weight  |  Tried   | Accepted | Acc. Ratio| Parameters" << std::endl;
        std::cout << "===============================================================================================================================" << std::endl;
        for (RbIterator<Move> it = moves.begin(); it != moves.end(); ++it)
        {
            it->printSummary(std::cout, current_period);
        }
        
        std::cout << std::endl;
        std::cout.flush();
    }
    
}


void Mcmc::replaceDag(const RbVector<Move> &mvs, const RbVector<Monitor> &mons)
{
    
    moves.clear();
    monitors.clear();
    
    // we need to replace the DAG nodes of the monitors and moves
    const std::vector<DagNode*>& model_nodes = model->getDagNodes();
    for (RbConstIterator<Move> it = mvs.begin(); it != mvs.end(); ++it)
    {
        
        Move *the_move = it->clone();
        std::vector<DagNode*> nodes = the_move->getDagNodes();
        for (std::vector<DagNode*>::const_iterator j = nodes.begin(); j != nodes.end(); ++j)
        {
            
            RevBayesCore::DagNode *the_node = *j;
            
            // error checking
            if ( the_node->getName() == "" )
            {
                throw RbException() << "Unable to connect move '" << the_move->getMoveName() << "' to DAG copy because variable name was lost"; 
            }
            
            DagNode* the_new_node = NULL;
            for (std::vector<DagNode*>::const_iterator k = model_nodes.begin(); k != model_nodes.end(); ++k)
            {
                if ( (*k)->getName() == the_node->getName() )
                {
                    the_new_node = *k;
                    break;
                }
            }
            // error checking
            if ( the_new_node == NULL )
            {
                throw RbException() << "Cannot find node with name '" << the_node->getName() << "' in the model but received a move working on it.";
            }
            
            // now swap the node
            the_move->swapNode( *j, the_new_node );
        }
        moves.push_back( *the_move );
        delete the_move;
    }
    
    for (RbConstIterator<Monitor> it = mons.begin(); it != mons.end(); ++it)
    {
        Monitor *the_monitor = it->clone();
        std::vector<DagNode*> nodes = the_monitor->getDagNodes();
        for (std::vector<DagNode*>::const_iterator j = nodes.begin(); j != nodes.end(); ++j)
        {
            
            RevBayesCore::DagNode *the_node = (*j);
            
            // error checking
            if ( the_node->getName() == "" )
            {
                throw RbException( "Unable to connect monitor to DAG copy because variable name was lost");
            }
            
            DagNode* the_new_node = NULL;
            for (std::vector<DagNode*>::const_iterator k = model_nodes.begin(); k != model_nodes.end(); ++k)
            {
                if ( (*k)->getName() == the_node->getName() )
                {
                    the_new_node = *k;
                    break;
                }
            }
            // error checking
            if ( the_new_node == NULL )
            {
                throw RbException() << "Cannot find node with name '" << the_node->getName() << "' in the model but received a monitor working on it.";
            }
            
            // now swap the node
            the_monitor->swapNode( *j, the_new_node );
        }
        monitors.push_back( *the_monitor );
        delete the_monitor;
        
    }
    
    resetVariableDagNodes();
}


void Mcmc::redrawStartingValues( void )
{
    
    std::vector<DagNode *> ordered_stoch_nodes = model->getOrderedStochasticNodes(  );
    for (std::vector<DagNode *>::iterator i=ordered_stoch_nodes.begin(); i!=ordered_stoch_nodes.end(); ++i)
    {
        DagNode *the_node = (*i);
        
        if ( the_node->isClamped() == false && the_node->isStochastic() == true )
        {

            the_node->redraw();
            the_node->reInitialized();
            
        }
        
        the_node->touch();
        
    }
    
    for (std::vector<DagNode *>::iterator i=ordered_stoch_nodes.begin(); i!=ordered_stoch_nodes.end(); ++i)
    {
        
        DagNode *the_node = (*i);
        the_node->keep();
        
    }
    
}


void Mcmc::removeMonitors( void )
{
    
    // just clear the vector
    monitors.clear();
    
}


/**
 * Reset the sampler.
 * We reset the counters of all moves.
 */
void Mcmc::reset( void )
{
    
    double moves_per_iteration = 0.0;
    for (RbIterator<Move> it = moves.begin(); it != moves.end(); ++it)
    {

        it->resetCounters();
        moves_per_iteration += it->getUpdateWeight();
        
    }

}


/**
 * Reset the currently monitored DAG nodes by extracting the DAG nodes from the StochasticVariable again
 * and store this in the set of DAG nodes.
 */
void Mcmc::resetVariableDagNodes( void )
{
    
    // for savety we empty our dag nodes
    variable_nodes.clear();
    
    if ( model != NULL )
    {
        // we only want to have each nodes once
        // this should by default happen by here we check again
        std::set<std::string> var_names;
        
        const std::vector<DagNode*> &n = model->getDagNodes();
        for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
        {
            
            DagNode *the_node = *it;
            
            // only non clamped variables
            if ( the_node->isClamped() == false )
            {
                if ( the_node->isStochastic() && the_node->isHidden() == false )
                {
                    const std::string &name = the_node->getName();
                    if ( var_names.find( name ) == var_names.end() )
                    {
                        variable_nodes.push_back( the_node );
                        var_names.insert( name );
                    }
                    else
                    {
#ifdef DEBUG_SEBASTIAN
                        std::cerr << "Trying to add variable with name '" << name << "' twice." << std::endl;
#endif
                    }
                    
                }
                
            }
            
        }
        
    }
    
}


/**
 * Set the active PID of this specific MCMC simulation.
 */
void Mcmc::setActivePIDSpecialized(size_t a, size_t n)
{
    
    // delegate the call to the model
    model->setActivePID(a,n);
    
    
    // tell each monitor
    for (size_t i=0; i < monitors.size(); ++i)
    {
        
        if ( process_active == true )
        {
            monitors[i].enable();
        }
        else
        {
            monitors[i].disable();
        }
        
    }
    
}


/**
 * Set if the current chain is the active chain.
 * Only active chains print to the monitors.
 */
void Mcmc::setChainActive(bool tf)
{
    chain_active = tf;
}


/**
 * Set the heat of the likelihood of the current chain.
 * This heat is used in posterior posterior MCMC algorithms to
 * heat the likelihood
 * The heat is passed to the moves for the accept-reject mechanism.
 */
void Mcmc::setChainLikelihoodHeat(double h)
{
    chain_likelihood_heat = h;
}


void Mcmc::setCheckpointFile(const path &f)
{
    checkpoint_file_name = f;
}


/**
 * Set the heat of the likelihood of the current chain.
 * This heat is used in power posterior MCMC algorithms to
 * heat the likelihood
 * The heat is passed to the moves for the accept-reject mechanism.
 */
void Mcmc::setLikelihoodHeat(double h)
{
    chain_likelihood_heat = h;
}


/**
 * Set the heat of the posterior of the current chain.
 * The heat of the posterior is used in the MC^3 algorithm.
 * The heat is passed to the moves for the accept-reject mechanism.
 */
void Mcmc::setChainPosteriorHeat(double h)
{
    chain_posterior_heat = h;
}


void Mcmc::setChainPriorHeat(double h)
{
    chain_prior_heat = h;
}


/**
 * Set the index of the current chain.
 */
void Mcmc::setChainIndex(size_t x)
{
    chain_idx = x;
}


/**
 * Set the model by delegating the model to the chains.
 */
void Mcmc::setModel( Model *m, bool redraw )
{
    // remember the old model
    Model * old_model = model;
    
    model = m;
    
    // we also need to replace the DAG nodes of our moves and monitors.
    RbVector<Move> tmp_moves = moves;
    RbVector<Monitor> tmp_monitors = monitors;
    replaceDag(tmp_moves, tmp_monitors);
    
    initializeMonitors();
    
    if ( redraw == true )
    {
        redrawStartingValues();
    }
    
    
    // free the old model
    delete old_model;
}


/**
 * Set the vector of moves of the current chain.
 */
void Mcmc::setMoves(const RbVector<Move> &mvs)
{
    moves = mvs;
}


void Mcmc::setMovesTuningInfo(const std::vector<tuningInfo> &mvs_ti)
{
    moves_tuningInfo = mvs_ti;
    
    if (moves_tuningInfo.size() != moves.size())
    {
        throw RbException( "The number of moves does not match the length of tuning information structures." );
    }
    
    for (size_t i = 0; i < moves.size(); ++i)
    {
        moves[i].setNumberTriedCurrentPeriod(moves_tuningInfo[i].num_tried_current_period);
        moves[i].setNumberTriedTotal(moves_tuningInfo[i].num_tried_total);
        moves[i].setNumberAcceptedCurrentPeriod(moves_tuningInfo[i].num_accepted_current_period);
        moves[i].setNumberAcceptedTotal(moves_tuningInfo[i].num_accepted_total);
        
        double tmp_tuningParameter = moves[i].getMoveTuningParameter();
        
        if ((std::isnan(moves_tuningInfo[i].tuning_parameter) == true && std::isnan(tmp_tuningParameter) == false) || (std::isnan(moves_tuningInfo[i].tuning_parameter) == false && std::isnan(tmp_tuningParameter) == true))
        {
            throw RbException( "The tunability of some moves changed." );
        }
        else if (std::isnan(tmp_tuningParameter) == false)
        {
            moves[i].setMoveTuningParameter(moves_tuningInfo[i].tuning_parameter);
        }
    }
    
}


void Mcmc::setScheduleType(const std::string &s)
{
    
    schedule_type = s;
}


/**
 * Start the monitors which will open the output streams.
 */
void Mcmc::startMonitors( size_t num_cycles, bool reopen )
{
    
    // Open the output file and print headers
    for (size_t i=0; i<monitors.size(); ++i)
    {
        
        // open filestream for each monitor
        monitors[i].openStream( reopen );
        
        // reset the monitor
        monitors[i].reset( num_cycles );
        
    }
    
}

/**
 * Write the header for each of the monitors.
 */
void Mcmc::writeMonitorHeaders( bool screen_monitor_only )
{
    
    // Open the output file and print headers
    for (size_t i=0; i<monitors.size(); ++i)
    {
        
        // if this chain is active, print the header
        if ( chain_active == true && process_active == true )
        {
            if ( monitors[i].isScreenMonitor() == true || screen_monitor_only == false )
            {
                monitors[i].printHeader();
            }
            
        }
        
    }
    
}


/**
 * Tune the sampler.
 * Here we just tune all the moves.
 */
void Mcmc::tune( void )
{
    
    // iterate over the moves
    for (RbIterator<Move> it=moves.begin(); it!=moves.end(); ++it)
    {
        // tune the move
        it->autoTune();
    }
    
}

