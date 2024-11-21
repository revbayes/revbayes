#include <cstddef>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "DagNode.h"
#include "HillClimber.h"
#include "MoveSchedule.h"
#include "RandomMoveSchedule.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RlUserInterface.h"
#include "SingleRandomMoveSchedule.h"
#include "SequentialMoveSchedule.h"
#include "MaximumLikelihoodEstimation.h"
#include "Model.h"
#include "Monitor.h"
#include "Move.h"
#include "RbConstIterator.h"
#include "RbConstIteratorImpl.h"
#include "RbIterator.h"
#include "RbIteratorImpl.h"
#include "RbVector.h"
#include "RbVectorImpl.h"


using namespace RevBayesCore;


/**
 * Constructor. We create an independent copy of the model and thus of all DAG nodes.
 * Someone might have wanted to run another HillClimber with different settings on the same model.
 * Thus we also create our own copies of the monitors and moves.
 *
 * \param[in]    m    The model containing all DAG nodes.
 * \param[in]    mvs  The vector of moves.
 * \param[in]    mons The vector of monitors.
 */
HillClimber::HillClimber(const Model& m, const RbVector<Move> &mvs, const RbVector<Monitor> &mons) : MaximumLikelihoodEstimation(),
    model( m.clone() ),
    monitors( mons ),
    moves( mvs ),
    schedule(NULL),
    scheduleType("random")
{
    // create an independent copy of the model, monitors and moves
    replaceDag(mvs,mons);

    initializeSampler();
    initializeMonitors();

}


/**
 * Copy constructor. For more details see the constructor.
 *
 * \param[in]    m    The MonteCarloSampler object to copy.
 */
HillClimber::HillClimber(const HillClimber &m) : MaximumLikelihoodEstimation(m),
    model( m.model->clone() ),
    monitors( m.monitors ),
    moves( m.moves ),
    schedule( NULL ),
    scheduleType( m.scheduleType )
{

    // temporary references
    const RbVector<Monitor>& mons = m.monitors;
    const RbVector<Move>& mvs = m.moves;


    // create an independent copy of the model, monitors and moves
    replaceDag(mvs,mons);

    initializeSampler();
    initializeMonitors();
}


/**
 * Destructor. Frees the DAG nodes (the model), moves, monitor and the move schedule.
 */
HillClimber::~HillClimber(void)
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
HillClimber& HillClimber::operator=(const HillClimber &m)
{
    MaximumLikelihoodEstimation::operator=( m );

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


///**
// * Add an extension to the name of the monitor.
// * We tell this to all our monitors.
// */
//void HillClimber::addFileMonitorExtension(const std::string &s, bool dir)
//{
//
//    // tell each monitor
//    for (RbIterator<Monitor> it=monitors.begin(); it!=monitors.end(); ++it)
//    {
//        it->addFileExtension( s, dir );
//    }
//
//}
//
//
//void HillClimber::addMonitor(const Monitor &m)
//{
//
//    monitors.push_back( m );
//
//}


/**
 * Disable all screen monitors. This means we simply delete it.
 */
void HillClimber::disableScreenMonitor( void )
{

    // tell each monitor
    for (size_t i=0; i < monitors.size(); ++i)
    {

        bool is = monitors[i].isScreenMonitor();
        if ( is == true )
        {
            monitors.erase( i );
            --i;
        }

    }

}


HillClimber* HillClimber::clone( void ) const
{

    return new HillClimber( *this );
}


/**
 * Finish the monitors which will close the output streams.
 */
void HillClimber::finishMonitors( void )
{

    // iterate over all monitors
    for (size_t i=0; i<monitors.size(); i++)
    {

        // if this chain is active, then close the stream
        if ( process_active == true )
        {
            monitors[i].closeStream();

        }

    }

}


/**
 * Get the model instance.
 */
Model& HillClimber::getModel( void )
{

    return *model;
}


/**
 * Get the model instance.
 */
const Model& HillClimber::getModel( void ) const
{

    return *model;
}


/**
 * Get the joint posterior probability of the current state for this model.
 * Note that the joint posterior is the true, unscaled and unheated value.
 */
double HillClimber::getModelLnProbability(bool likelihood_only)
{
    const std::vector<DagNode*> &n = model->getDagNodes();
    double pp = 0.0;
    for (std::vector<DagNode*>::const_iterator it = n.begin(); it != n.end(); ++it)
    {
        DagNode *the_node = *it;
        if ( likelihood_only == true && the_node->isClamped() == true )
        {
            pp += the_node->getLnProbability();
        }
        else
        {
          pp += the_node->getLnProbability();
        }

    }
    return pp;
}


/**
 * Get the vector of monitors for this sampler.
 */
RbVector<Monitor>& HillClimber::getMonitors(void)
{
    return monitors;
}


/**
 * Get the vector of moves for this sampler.
 */
RbVector<Move>& HillClimber::getMoves(void)
{
    return moves;
}


/**
 * Get a const-reference move-schedule for this sampler.
 */
const MoveSchedule& HillClimber::getSchedule(void) const
{
    return *schedule;
}


/**
 * Get a non-const reference to the move-schedule of this sampler.
 */
MoveSchedule& HillClimber::getSchedule(void)
{
    return *schedule;
}

/**
 * Get the schedule type of this sampler.
 */
const std::string& HillClimber::getScheduleType( void ) const
{
    return scheduleType;
}


std::string HillClimber::getStrategyDescription( void ) const
{

    std::string description = "";
    std::stringstream stream;
    if ( scheduleType == "single" )
    {
        stream << "The simulator uses " << moves.size() << " different moves, with a single move picked randomly per iteration" << std::endl;
    }
    else if ( scheduleType == "random" )
    {
        stream << "The simulator uses " << moves.size() << " different moves in a random move schedule with " << schedule->getNumberMovesPerIteration() << " moves per iteration" << std::endl;
    }
    else if ( scheduleType == "sequential" )
    {
        stream << "The simulator uses " << moves.size() << " different moves in a sequential move schedule with " << schedule->getNumberMovesPerIteration() << " moves per iteration" << std::endl;
    }
    description = stream.str();

    return description;
}


bool HillClimber::hasConverged(double min_acceptance_ratio)
{
    bool converged = true;

    for (size_t i=0; i<moves.size() && converged == true; ++i)
    {
        size_t num_tried    = moves[i].getNumberTriedCurrentPeriod();
        size_t num_accepted = moves[i].getNumberAcceptedCurrentPeriod();

        if ( num_tried > 0)
        {
            bool tmp = ( min_acceptance_ratio > (double(num_accepted)/double(num_tried)) );
            converged &= tmp;
        }

    }

    return converged;
}


void HillClimber::initializeSampler( void )
{

    std::vector<DagNode *>& dagNodes = model->getDagNodes();
    std::vector<DagNode *> orderedStochNodes = model->getOrderedStochasticNodes(  );

    // Get rid of previous move schedule, if any
    if ( schedule )
    {
        delete schedule;
    }
    schedule = NULL;

    // Get initial lnProbability of model

    // first we touch all nodes so that the likelihood is dirty
    for (std::vector<DagNode *>::iterator i=dagNodes.begin(); i!=dagNodes.end(); i++)
    {

        DagNode *the_node = *i;
        the_node->setMcmcMode( true );
        the_node->setPriorOnly( false );
        the_node->touch();

    }


    int numTries    = 0;
    int maxNumTries = 100;
    double lnProbability = 0.0;
    for ( ; numTries < maxNumTries; numTries ++ )
    {
        // a flag if we failed to find a valid starting value
        bool failed = false;

        lnProbability = 0.0;
        for (std::vector<DagNode *>::iterator i=dagNodes.begin(); i!=dagNodes.end(); i++)
        {
            DagNode* node = (*i);
            node->touch();

            double lnProb = node->getLnProbability();

            if ( !RbMath::isAComputableNumber(lnProb) )
            {
                std::stringstream ss;
                ss << "Could not compute lnProb for node " << node->getName() << "." << std::endl;
                node->printValue( ss, "," );
                ss << std::endl;
                RBOUT( ss.str() );

                // set the flag
                failed = true;

                break;
            }
            lnProbability += lnProb;

        }

        // now we keep all nodes so that the likelihood is stored
        for (std::vector<DagNode *>::iterator i=dagNodes.begin(); i!=dagNodes.end(); i++)
        {
            (*i)->keep();
        }

        if ( failed == true )
        {
            RBOUT( "Drawing new initial states ... " );
            for (std::vector<DagNode *>::iterator i=orderedStochNodes.begin(); i!=orderedStochNodes.end(); i++)
            {

                if ( !(*i)->isClamped() && (*i)->isStochastic() )
                {
                    (*i)->redraw();
                    (*i)->reInitialized();

                }
                else if ( (*i)->isClamped() )
                {
                    // make sure that the clamped node also recompute their probabilities
                    (*i)->reInitialized();
                    (*i)->touch();
                }

            }
        }
        else
        {
            break;
        }

    }

    if ( numTries == maxNumTries )
    {
        std::stringstream msg;
        msg << "Unable to find a starting state with computable probability";
        if ( numTries > 1 )
        {
            msg << " after " << numTries << " tries";
        }
        throw RbException( msg.str() );

    }

    // Create the move scheduler
    if ( scheduleType == "sequential" )
    {
        schedule = new SequentialMoveSchedule( &moves );
    }
    else if ( scheduleType == "single" )
    {
        schedule = new SingleRandomMoveSchedule( &moves );
    }
    else
    {
        schedule = new RandomMoveSchedule( &moves );
    }

    generation = 0;
}



void HillClimber::initializeMonitors(void)
{
    for (size_t i=0; i<monitors.size(); i++)
    {
        monitors[i].setModel( model );
    }
}


void HillClimber::monitor(unsigned long g)
{

    if ( process_active == true )
    {
        // Monitor
        for (size_t i = 0; i < monitors.size(); i++)
        {
            monitors[i].monitor( g );
        }
    }

}


void HillClimber::nextCycle( void )
{

    size_t proposals = size_t( round( schedule->getNumberMovesPerIteration() ) );
    for (size_t i=0; i<proposals; i++)
    {

        // Get the move
        Move& the_move = schedule->nextMove( generation );

        // Perform the move
        the_move.performHillClimbingStep( 1.0, 1.0);

    }


    // advance gen cycle
    generation++;

}



void HillClimber::replaceDag(const RbVector<Move> &mvs, const RbVector<Monitor> &mons)
{

    moves.clear();
    monitors.clear();

    // we need to replace the DAG nodes of the monitors and moves
    const std::vector<DagNode*>& modelNodes = model->getDagNodes();
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
                throw RbException( "Unable to connect move '" + the_move->getMoveName() + "' to DAG copy because variable name was lost");
            }

            DagNode* theNewNode = NULL;
            for (std::vector<DagNode*>::const_iterator k = modelNodes.begin(); k != modelNodes.end(); ++k)
            {
                if ( (*k)->getName() == the_node->getName() )
                {
                    theNewNode = *k;
                    break;
                }
            }
            // error checking
            if ( theNewNode == NULL )
            {
                throw RbException() << "Cannot find node with name '" << the_node->getName() << "' in the model but received a move working on it.";
            }

            // now swap the node
            the_move->swapNode( *j, theNewNode );
        }
        moves.push_back( *the_move );
        delete the_move;
    }

    for (RbConstIterator<Monitor> it = mons.begin(); it != mons.end(); ++it)
    {
        Monitor *theMonitor = it->clone();
        std::vector<DagNode*> nodes = theMonitor->getDagNodes();
        for (std::vector<DagNode*>::const_iterator j = nodes.begin(); j != nodes.end(); ++j)
        {

            RevBayesCore::DagNode *the_node = (*j);

            // error checking
            if ( the_node->getName() == "" )
            {
                throw RbException( "Unable to connect monitor to DAG copy because variable name was lost");
            }

            DagNode* theNewNode = NULL;
            for (std::vector<DagNode*>::const_iterator k = modelNodes.begin(); k != modelNodes.end(); ++k)
            {
                if ( (*k)->getName() == the_node->getName() )
                {
                    theNewNode = *k;
                    break;
                }
            }
            // error checking
            if ( theNewNode == NULL )
            {
                throw RbException() << "Cannot find node with name '" << the_node->getName() << "' in the model but received a monitor working on it.";
            }

            // now swap the node
            theMonitor->swapNode( *j, theNewNode );
        }
        monitors.push_back( *theMonitor );
        delete theMonitor;

    }

}


void HillClimber::removeMonitors( void )
{

    // just clear the vector
    monitors.clear();

}


/**
 * Reset the sampler.
 * We reset the counters of all moves.
 */
void HillClimber::reset( void )
{

    double movesPerIteration = 0.0;
    for (RbIterator<Move> it = moves.begin(); it != moves.end(); ++it)
    {

        it->resetCounters();
        movesPerIteration += it->getUpdateWeight();

    }

}


/**
 * Set the active PID of this specific HillClimber simulation.
 */
void HillClimber::setActivePIDSpecialized(size_t a, size_t n)
{

    // delegate the call to the model
    model->setActivePID(a,n);
}


/**
 * Set the model by delegating the model to the chains.
 */
void HillClimber::setModel( Model *m )
{

    model = m;

    // we need to set the new model in the monitors as well
    // which is done in initializeMonitors
    initializeMonitors();

}


void HillClimber::setScheduleType(const std::string &s)
{

    scheduleType = s;
}


/**
 * Start the monitors which will open the output streams.
 */
void HillClimber::startMonitors( size_t num_cycles, bool reopen )
{

    // Open the output file and print headers
    for (size_t i=0; i<monitors.size(); ++i)
    {

        // reset the monitor
        monitors[i].reset( num_cycles );

        // if this chain is active, print the header
        if ( process_active == true )
        {
            monitors[i].openStream( reopen );
        }

    }

}


/**
 * Tune the sampler.
 * Here we just tune all the moves.
 */
void HillClimber::tune( void )
{

    // iterate over the moves
    for (RbIterator<Move> it=moves.begin(); it!=moves.end(); ++it)
    {
        // tune the move
        it->autoTune();
    }

}


/**
 * Start the monitors which will open the output streams.
 */
void HillClimber::writeMonitorHeaders( void )
{

    // Open the output file and print headers
    for (size_t i=0; i<monitors.size(); ++i)
    {

        // if this chain is active, print the header
        if ( process_active == true )
        {

            monitors[i].printHeader();

        }

    }

}
