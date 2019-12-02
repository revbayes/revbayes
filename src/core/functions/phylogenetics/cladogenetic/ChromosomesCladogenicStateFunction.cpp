#include "ChromosomesCladogenicStateFunction.h"
#include "MatrixReal.h"
#include "RbException.h"


using namespace RevBayesCore;


ChromosomesCladogenicStateFunction::ChromosomesCladogenicStateFunction(const TypedDagNode< Simplex > *ep, unsigned mc):
    TypedFunction<CladogeneticProbabilityMatrix>( new CladogeneticProbabilityMatrix( mc + 1 ) ),
    eventProbs( ep ),
    maxChromo(mc),
    numEventTypes( (unsigned)ep->getValue().size() )
{
    addParameter( eventProbs );
    
    // since the transition probability matrix will be very large (maxChromo by maxChromo^2) but
    // only sparsely filled, the PhyloCTMCClado uses an event map instead
    // of the full matrix for efficiency
    buildEventMap();
    update();
}


ChromosomesCladogenicStateFunction::~ChromosomesCladogenicStateFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


void ChromosomesCladogenicStateFunction::buildEventMap( void ) {
    
    eventMapCounts.resize(maxChromo + 1, std::vector<unsigned>(numEventTypes, 0));

    // for each ancestor state build a map of the possible events in the structure:
    // pair< [ancestor_state, daughter_1_state, daughter_2_state], transition_probability >
    std::vector<unsigned> idx(3);
    for (unsigned i = 1; i <= maxChromo; i++)
    {
        // set ancestor state
        idx[0] = i;

        // loop through all possible events
        for (unsigned j = 0; j < numEventTypes; j++)
        {
            // if both daughters have same number chromsomes as ancestor
            if (j == NO_CHANGE)
            {
                idx[1] = i;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(NO_CHANGE));
                eventMapCounts[ i ][ NO_CHANGE ] += 1;
                eventMapProbs[ idx ] = 0.0;
            }
            // if one of the daughters gains a chromosome
            else if (j == FISSION && i + 1 <= maxChromo)
            {
                idx[1] = i + 1;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(FISSION));
                eventMapCounts[ i ][ FISSION ] += 1;
                eventMapProbs[ idx ] = 0.0;

                idx[1] = i;
                idx[2] = i + 1;
                eventMapTypes[ idx ].push_back(unsigned(FISSION));
                eventMapCounts[ i ][ FISSION ] += 1;
                eventMapProbs[ idx ] = 0.0;
            }
            // if one of the daughters loses a chromosome
            else if ( j == FUSION && i > 1 )
            {
                idx[1] = i - 1;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(FUSION));
                eventMapCounts[ i ][ FUSION ] += 1;
                eventMapProbs[ idx ] = 0.0;

                idx[1] = i;
                idx[2] = i - 1;
                eventMapTypes[ idx ].push_back(unsigned(FUSION));
                eventMapCounts[ i ][ FUSION ] += 1;
                eventMapProbs[ idx ] = 0.0;
            }
            // if one of the daughters undergoes polyploidizations
            else if ( j == POLYPLOIDIZATION && i * 2 <= maxChromo )
            {
                idx[1] = i * 2;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(POLYPLOIDIZATION));
                eventMapCounts[ i ][ POLYPLOIDIZATION ] += 1;
                eventMapProbs[ idx ] = 0.0;

                idx[1] = i;
                idx[2] = i * 2;
                eventMapTypes[ idx ].push_back(unsigned(POLYPLOIDIZATION));
                eventMapCounts[ i ][ POLYPLOIDIZATION ] += 1;
                eventMapProbs[ idx ] = 0.0;
            }
            // if one of the daughters undergoes demipolyploidizations
            else if ( j == DEMIPOLYPLOIDIZATION && i * 1.5 <= maxChromo && i > 1 )
            {
                if (i % 2 == 0)
                {
                    idx[1] = i * 1.5;
                    idx[2] = i;
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMapProbs[ idx ] = 0.0;

                    idx[1] = i;
                    idx[2] = i * 1.5;
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMapProbs[ idx ] = 0.0;
                }
                else
                {
                    // round down
                    idx[1] = (unsigned)( (double)i * 1.5 - 0.5 );
                    idx[2] = i;
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMapProbs[ idx ] = 0.0;

                    idx[1] = i;
                    idx[2] = (unsigned)( (double)i * 1.5 - 0.5 );
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMapProbs[ idx ] = 0.0;
                    
                    if ( i * 1.5 + 0.5 <= maxChromo )
                    {
                        // round up
                        idx[1] = i;
                        idx[2] = (unsigned)( (double)i * 1.5 + 0.5 );
                        eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                        eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                        eventMapProbs[ idx ] = 0.0;
                        
                        idx[1] = (unsigned)( (double)i * 1.5 + 0.5 );
                        idx[2] = i;
                        eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                        eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1; 
                        eventMapProbs[ idx ] = 0.0;
                    }
                }
            }
        }
    }
}


ChromosomesCladogenicStateFunction* ChromosomesCladogenicStateFunction::clone( void ) const
{
    return new ChromosomesCladogenicStateFunction( *this );
}

double ChromosomesCladogenicStateFunction::computeDataAugmentedCladogeneticLnProbability(const std::vector<BranchHistory*>& histories,
                                                                                              size_t node_index,
                                                                                              size_t left_index,
                                                                                              size_t right_index ) const
{
    throw RbException("ChromosomesCladogenicStateFunction::computeDataAugmentedCladogeneticLnProbability is not currently implemented.");
    double lnP = 0.0;
    return lnP;
    
}


std::map< std::vector<unsigned>, double >  ChromosomesCladogenicStateFunction::getEventMap(double t)
{
    return eventMapProbs;
}

const std::map< std::vector<unsigned>, double >&  ChromosomesCladogenicStateFunction::getEventMap(double t) const
{
    return eventMapProbs;
}


void ChromosomesCladogenicStateFunction::update( void )
{
    // reset the transition matrix
    delete value;
    value = new CladogeneticProbabilityMatrix( maxChromo + 1 );

    const std::vector<double>& ep = eventProbs->getValue();
    
    // normalize the transition probabilities
    // for each row of the transition matrix
    for (unsigned i = 1; i <= maxChromo; i++)
    {
        // loop through all the mapped events for this start state
        // and get the sum of probabilities
        double prob_sum = 0.0;
        std::map<std::vector<unsigned>, std::vector<unsigned> >::iterator it;
        for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
        {
            const std::vector<unsigned>& idx = it->first;
            if (idx[0] == i)
            {
                double event_prob = 0.0;
                std::vector<unsigned> event_types = it->second;
                for (size_t e = 0; e < event_types.size(); e++)
                {
                    // reset all the probs to 0.0
                    eventMapProbs[ idx ] = 0.0;
//                    (*value)[ idx[0] ][ (maxChromo + 1) * idx[1] + idx[2] ] = 0.0;
                    
                    // check for NaN values
                    if (ep[ event_types[e] ] == ep[ event_types[e] ])
                    {
                        event_prob = ep[ event_types[e] ];
                    }

                    // normalize for all possible instances of this event type
                    prob_sum += event_prob / eventMapCounts[ i ][ event_types[e] ];
                }
            } 
        }
        // now normalize the probabilities so they sum to 1.0
        for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
        {
            const std::vector<unsigned>& idx = it->first;
            if (idx[0] == i)
            {
                std::vector<unsigned> event_types = it->second;
                for (size_t e = 0; e < event_types.size(); e++)
                {
                    double event_prob = 0.0;
                    
                    // check for NaN values
                    if (ep[ event_types[e] ] == ep[ event_types[e] ])
                    {
                        event_prob = ep[ event_types[e] ];
                    }

                    // normalize for all possible instances of this event type
                    double v = ( event_prob / eventMapCounts[ i ][ event_types[e] ] ) / prob_sum;

                    // save the probability in the transition matrix
//                    (*value)[ idx[0] ][ (maxChromo + 1) * idx[1] + idx[2] ] += v;

                    // save the probability in the event map
                    eventMapProbs[ idx ] += v;
                }
            } 
        }
    }
    value->setEventMap(eventMapProbs);
}


void ChromosomesCladogenicStateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == eventProbs)
    {
        eventProbs = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
}

