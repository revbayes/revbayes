#include "StartingTreeSimulator.h"

#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <iosfwd>
#include <iterator>
#include <map>
#include <string>
#include <type_traits>
#include <utility>

#include "DistributionExponential.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"
#include "TopologyNode.h"
#include "Clade.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "Tree.h"

using namespace RevBayesCore;


StartingTreeSimulator::StartingTreeSimulator( void )
{
    
}


/**
 *
 */
Tree* StartingTreeSimulator::simulateTree( const std::vector<Taxon> &taxa, const std::vector<Clade> &constraints ) const
{
    std::vector<Clade> my_constraints = constraints;
    
    // the time tree object (topology & times)
    Tree *psi = new Tree();
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( true );
    
    size_t num_taxa = taxa.size();
    
    // add a clounter variable of how many missing taxa we have already added
    size_t n_added_missing_taxa = 0;
    
    // create the tip nodes
    std::vector<TopologyNode*> nodes;
    for (size_t i=0; i<num_taxa; ++i)
    {
        
        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa[i], i );
        
        // set the age of this tip node
        node->setAge( taxa[i].getAge() );
        
        // add the new node to the list
        nodes.push_back( node );
        
    }
    
    // we need a sorted vector of constraints, starting with the smallest
    std::vector<Clade> sorted_clades;
    
    for (size_t i = 0; i < my_constraints.size(); ++i)
    {
        
        // set the ages of each of the taxa in the constraint
        for (size_t j = 0; j < my_constraints[i].size(); ++j)
        {
            for (size_t k = 0; k < num_taxa; ++k)
            {
                if ( taxa[k].getName() == my_constraints[i].getTaxonName(j) )
                {
                    my_constraints[i].setTaxonAge(j, taxa[k].getAge());
                    break;
                }
            }
        }
        
        // set ages for optional constraints
        if ( my_constraints[i].hasOptionalConstraints() == true )
        {
            std::vector<Clade> optional_constraints = my_constraints[i].getOptionalConstraints();
            for (size_t k = 0; k < optional_constraints.size(); k++)
            {
                for (size_t opt_taxon_idx = 0; opt_taxon_idx < optional_constraints[k].size(); opt_taxon_idx++)
                {
                    for (size_t full_taxon_idx = 0; full_taxon_idx < num_taxa; full_taxon_idx++)
                    {
                        if ( taxa[full_taxon_idx].getName() == optional_constraints[k].getTaxonName(opt_taxon_idx) )
                        {
                        
                            optional_constraints[k].setTaxonAge(opt_taxon_idx, taxa[full_taxon_idx].getAge());
                            break;
                        }
                    }
                }
            
            }
        
            my_constraints[i].setOptionalConstraints( optional_constraints );
        } // finished with the optional constrains, if they exist
        
        // populate sorted clades vector
        if ( my_constraints[i].size() > 1 && my_constraints[i].size() <= num_taxa )
        {
            
            if ( my_constraints[i].hasOptionalConstraints() == true )
            {
                std::vector<Clade> optional_constraints = my_constraints[i].getOptionalConstraints();
                size_t idx = (size_t)( GLOBAL_RNG->uniform01() * optional_constraints.size() );
                sorted_clades.push_back( optional_constraints[idx] );
            }
            else
            {
                sorted_clades.push_back( my_constraints[i] );
            }
        }
        
    }

    // DO WE NEED TO SORT THE TAXA?
    std::sort(sorted_clades.begin(), sorted_clades.end(), cladeBefore);

    if(sorted_clades.back().size() < num_taxa) {
        // create a clade that contains all species
        Clade all_species = Clade(taxa);
        all_species.setAge( RbConstants::Double::inf );
        sorted_clades.push_back(all_species);
    }
    
    std::vector<Clade> virtual_taxa;
    int i = -1;
    for (std::vector<Clade>::iterator it = sorted_clades.begin(); it != sorted_clades.end(); it++)
    {
        // ignore negative clade constraints during simulation
        if ( it->isNegativeConstraint() == true )
        {
            continue;
        }
        
        ++i;
        const Clade &c = *it;
        std::vector<Taxon> taxa = c.getTaxa();
        std::vector<Clade> clades;
        
        int j = i;
        std::vector<Clade>::reverse_iterator jt(it);
        for (; jt != sorted_clades.rend(); jt++)
        {
            // ignore negative clade constraints during simulation
            if ( jt->isNegativeConstraint() == true )
            {
                continue;
            }
            
            j--;
            const Clade &c_nested = *jt;
            std::vector<Taxon> taxa_nested = c_nested.getTaxa();
            
            bool found_all = true;
            bool found_some = false;
            for (size_t k = 0; k < taxa_nested.size(); ++k)
            {
                std::vector<Taxon>::iterator kt = std::find(taxa.begin(), taxa.end(), taxa_nested[k]);
                if ( kt != taxa.end() )
                {
                    taxa.erase( kt );
                    found_some = true;
                }
                else
                {
                    found_all = false;
                }
                
            }
            
            if ( found_all == true )
            {
                //                c.addTaxon( virtual_taxa[j] );
                //                taxa.push_back( virtual_taxa[j] );
                clades.push_back( virtual_taxa[j] );
            }
            
            if ( found_all == false && found_some == true )
            {
                throw RbException("Cannot simulate tree: conflicting monophyletic clade constraints. Check that all clade constraints are properly nested.");
            }
            
        }
        
        
        std::set<TopologyNode*> nodes_in_clade;
        
        
        for (size_t k = 0; k < taxa.size(); ++k)
        {
            Clade tmp_clade = Clade( taxa[k] );
            tmp_clade.setAge( taxa[k].getAge() );
            clades.push_back( tmp_clade );
        }
        
        for (size_t k = 0; k < clades.size(); ++k)
        {
            for (size_t j = 0; j < nodes.size(); ++j)
            {
                if (nodes[j]->getClade() == clades[k])
                {
                    nodes_in_clade.insert( nodes[j] );
                    nodes.erase( nodes.begin()+j );
                    break;
                }
                
            }
            
        }
        
        
        // here we need to start adding our "fake" tips
        for ( int index_missing_species = 0; index_missing_species < c.getNumberMissingTaxa(); ++index_missing_species)
        {
            ++n_added_missing_taxa;
            TopologyNode* missing_node = new TopologyNode("Missing_Taxon_" + StringUtilities::to_string(n_added_missing_taxa) );
            missing_node->setAge( 0.0 );
            nodes_in_clade.insert( missing_node );
        }
        
        simulateClade(nodes_in_clade, c.getAge());
        TopologyNode *remaining_node = *(nodes_in_clade.begin());
        nodes.push_back( remaining_node );
        
        std::vector<Taxon> v_taxa;
        remaining_node->getTaxa(v_taxa);
        Clade new_clade = Clade(v_taxa);
        remaining_node -> setAge(c.getAge());
        new_clade.setAge( remaining_node->getAge() );
        virtual_taxa.push_back( new_clade );
    }
    
    TopologyNode *root = nodes[0];
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    return psi;
}


void StartingTreeSimulator::simulateClade( std::set<TopologyNode*> &nodes, double max_clade_age) const
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    size_t num_taxa = nodes.size();
    
    // retrieve the times of any serially sampled tips
    std::vector<double> serial_times;
    size_t num_taxa_at_present = 0;
    
    double min_age = RbConstants::Double::inf;
    for (std::set<TopologyNode*>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        TopologyNode *this_node = *it;
        double a = this_node->getAge();
        if ( min_age > a )
        {
            min_age = a;
        }
    }

    double max_tip_age = 0;
    std::set<TopologyNode*> current_nodes;
    std::multimap<double,TopologyNode*> serial_nodes;
    for (std::set<TopologyNode*>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        TopologyNode *this_node = *it;
        double a = this_node->getAge();
        if ( this_node->isTip() == true )
        {
            a = this_node->getTaxon().getAgeRange().getMax();
        }
        if ( a > min_age )
        {
            serial_times.push_back(a);
            serial_nodes.insert( std::pair<double, TopologyNode*>(a, this_node) );
            if(a > max_tip_age) max_tip_age = a;
        }
        else
        {
            ++num_taxa_at_present;
            current_nodes.insert( this_node );
        }
    }

    if(max_tip_age >= max_clade_age) throw RbException() << "Specified clade age " << max_clade_age << " lower than max tip age";
    
    size_t index_serial_time = 0;
    if (num_taxa_at_present == num_taxa)
    {
        serial_times.push_back(max_clade_age);
    }
    else
    {
        std::sort(serial_times.begin(), serial_times.end());
        serial_times.push_back(max_clade_age);
    }
    
    // now simulate the ages
        
    // j is the number of active lineages at the current time
    size_t j = num_taxa_at_present;
    //    double theta = Ne->getValue();
    double theta = 1.0;
    
    // the current age of the process
    double sim_age = min_age;
    
    // draw a time for each speciation event condition on the time of the process
    size_t num_events_to_draw = nodes.size() - 1;
    while ( num_events_to_draw > 0)
    {
        bool valid = false;
        do
        {
            double nPairs = j * (j-1) / 2.0;
            double lambda = nPairs / theta;
            double u = RbStatistics::Exponential::rv( lambda, *rng);

            if(sim_age + u > max_clade_age && j > 1) {
                valid = false;
                theta /= 2; //we're simulating branch lengths that are too std::int64_t, so decreasing Ne
            } else {
                sim_age += u;
                valid = sim_age < serial_times[index_serial_time] && j > 1;
                if ( valid == false )
                {
                    // If j is 1 and we are still simulating coalescent events, we have >= 1 serial sample left to coalesce.
                    // There are no samples to coalesce now, but we cannot exit, thus, we advance to the next serial sample
                    // Alternately, when we cross a serial sampling time, the number of active lineages changes
                    // it is necessary to discard any "excess" time, which is drawn from an incorrect distribution
                    // then we can draw a new time according to the correct number of active lineages.
                    // Either we advance or go back, but in both situations we set the time to the current serial sample.
                    sim_age = serial_times[index_serial_time];
                    ++index_serial_time;
                    ++j;
                    
                    std::multimap<double,TopologyNode*>::iterator it = serial_nodes.begin();
                    current_nodes.insert( it->second );
                    serial_nodes.erase( it );
                    
                }
            }
        } while ( valid == false );
        
        size_t left_index = size_t( current_nodes.size() * rng->uniform01() );
        std::set<TopologyNode*>::iterator left_it = current_nodes.begin();
        std::advance(left_it, left_index);
        TopologyNode *left = *left_it;
        current_nodes.erase( left_it );
        
        size_t right_index = size_t( current_nodes.size() * rng->uniform01() );
        std::set<TopologyNode*>::iterator right_it = current_nodes.begin();
        std::advance(right_it, right_index);
        TopologyNode *right = *right_it;
        current_nodes.erase( right_it );
        
        TopologyNode *parent = new TopologyNode();
        parent->addChild( left );
        parent->addChild( right );
        left->setParent( parent );
        right->setParent( parent );
        
        parent->setAge( sim_age );
        current_nodes.insert( parent );

        --j;
        --num_events_to_draw;
        
    }
    
    nodes.clear();
    nodes.insert( *(current_nodes.begin()) );
}
