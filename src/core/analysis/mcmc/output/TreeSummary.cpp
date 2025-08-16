#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <optional>
#include <range/v3/all.hpp>

#include "NewickConverter.h"
#include "ProgressBar.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbVectorUtilities.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "TopologyNode.h"
#include "TreeSummary.h"
#include "TreeUtilities.h"
#include "Clade.h"
#include "RbBitSet.h"
#include "Taxon.h"
#include "TaxonMap.h"
#include "Trace.h"
#include "TraceTree.h"
#include "Tree.h"

using namespace RevBayesCore;

namespace views = ranges::views;

/*
 * TreeSummary constructor
 */
TreeSummary::TreeSummary( TraceTree* t, bool c ) :
    clock( c ),
    rooted( c )
{
    traces.push_back(t);
}


/*
 * TreeSummary constructor
 */
TreeSummary::TreeSummary( std::vector<TraceTree* > t, bool c ) :
    traces(t),
    clock( c ),
    rooted( c )
{
    if( traces.empty() )
    {
        throw(RbException("Tree summary requires at least one tree trace"));
    }

    std::vector<std::string> tip_names = traces.front()->objectAt(0).getTipNames();
    std::sort(tip_names.begin(),tip_names.end());

    for(auto& trace: traces)
    {
        std::vector<std::string> t = trace->objectAt(0).getTipNames();
        std::sort(t.begin(),t.end());

        if( t != tip_names )
        {
            throw(RbException("Cannot summarize tree traces with mismatched tip names"));
        }
    }
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
TreeSummary* TreeSummary::clone(void) const
{

    return new TreeSummary(*this);
}


void TreeSummary::annotateTree( Tree &tree, AnnotationReport report, bool verbose )
{
    summarize( verbose );

    RBOUT("Annotating tree ...");

    std::string newick;

    if ( report.conditional_tree_ages )
    {
        Tree* tmp_tree = NULL;
        if ( clock == true )
        {
            tmp_tree = TreeUtilities::convertTree( tree );
        }
        else
        {
            tmp_tree = tree.clone();
        }

        if ( tmp_tree->isRooted() == false && rooted == false )
        {
            if ( outgroup )
            {
                tmp_tree->reroot( *outgroup, false, true );
            }
            else
            {
                std::vector<std::string> tip_names = traces.front()->objectAt(0).getTipNames();
                std::sort(tip_names.begin(),tip_names.end());
                std::string outgrp = tip_names[0];
                tmp_tree->reroot( outgrp, false, true );
            }
        }
        else if ( tmp_tree->isRooted() != rooted )
        {
            throw(RbException("Rooting of input tree differs from the tree sample"));
        }

        newick = tmp_tree->getPlainNewickRepresentation();

        delete tmp_tree;

        if ( tree_clade_ages.find(newick) == tree_clade_ages.end() )
        {
            throw(RbException("Could not find input tree in tree sample"));
        }
    }

    const std::vector<TopologyNode*> &nodes = tree.getNodes();

    double totalSamples = sampleSize(true);

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        TopologyNode* n = nodes[i];

        Clade clade = n->getClade();
        Split split( clade.getBitRepresentation(), clade.getMrca(), rooted);

        // annotate clade posterior prob
        if ( ( !n->isTip() || ( n->isRoot() && !clade.getMrca().empty() ) ) && report.clade_probs )
        {
            double pp = cladeProbability( clade, false );
            n->addNodeParameter("posterior",pp);
        }

        // are we using sampled ancestors?
        if ( sampled_ancestor_counts.empty() == false )
        {
            double saFreq = sampled_ancestor_counts[n->getTaxon()];

            // annotate sampled ancestor prob
            if ( ((n->isTip() && n->isFossil()) || saFreq > 0) && report.sampled_ancestor_probs )
            {
                n->addNodeParameter("sampled_ancestor", saFreq / totalSamples);
            }
        }

        // annotate conditional clade probs and get node ages
        std::vector<double> node_ages;

        if ( !n->isRoot() )
        {
            Clade parent_clade = n->getParent().getClade();
            Split parent_split = Split( parent_clade.getBitRepresentation(), parent_clade.getMrca(), rooted);

            std::map<Split, std::vector<double> >& condCladeAges = conditional_clade_ages[parent_split];
            node_ages = report.conditional_clade_ages ? condCladeAges[split] : clade_ages[split];

            // annotate CCPs
            if ( !n->isTip() && report.conditional_clade_probs )
            {
                double parentCladeFreq = splitCount( parent_split );
                double ccp = condCladeAges[split].size() / parentCladeFreq;
                n->addNodeParameter("ccp",ccp);
            }
        }
        else
        {
            node_ages = clade_ages[split];
        }

        if ( report.conditional_tree_ages )
        {
            node_ages = tree_clade_ages[newick][split];
        }

        // set the node ages/branch lengths
        if ( report.node_ages )
        {
            double age = 0.0;

            if ( report.mean_node_ages )
            {
                // finally, we compute the mean conditional age
                for (size_t i = 0; i<node_ages.size(); ++i)
                {
                    age += node_ages[i];
                }
                age /= node_ages.size();
            }
            else // use median
            {

                size_t idx = node_ages.size() / 2;
                std::sort( node_ages.begin(), node_ages.end() );
                if (node_ages.size() % 2 == 1)
                {
                    age = node_ages[idx];
                }
                else
                {
                    age = (node_ages[idx-1] + node_ages[idx]) / 2;
                }

            }

            // finally, we set the age/length
            if ( clock )
            {
                n->setAge( age );
            }
            else
            {
                n->setBranchLength( age );
            }
        }

        // annotate the HPD node age intervals
        if ( report.node_ages_HPD )
        {
            //node_ages = cladeAges[c];

            std::sort(node_ages.begin(), node_ages.end());

            size_t total_branch_lengths = node_ages.size();
            double min_range = std::numeric_limits<double>::max();

            size_t interval_start = 0;
            double lower = node_ages[(int)(0.5 * (double)total_branch_lengths)];
            double upper = node_ages[(int)(0.5 * (double)total_branch_lengths)];

            int interval_size = (int)(report.node_ages_HPD * (double)total_branch_lengths);
            // we need to make sure that we sampled more than one age
            if ( interval_size > 1 )
            {

                // find the smallest interval that contains x% of the samples
                for (size_t j = 0; j <= (total_branch_lengths - interval_size); j++)
                {
                    double temp_lower = node_ages[j];
                    double temp_upper = node_ages[j + interval_size - 1];
                    double temp_range = std::fabs(temp_upper - temp_lower);
                    if (temp_range < min_range)
                    {
                        min_range = temp_range;
                        interval_start = j;
                    }

                }
                lower = node_ages[interval_start];
                upper = node_ages[interval_start + interval_size - 1];

            }

            // make node age annotation
            std::string interval = "{" + StringUtilities::toString(lower)
            + "," + StringUtilities::toString(upper) + "}";

            if ( clock == true )
            {
                if ( n->isTip() == false || ( ( n->isFossil() || upper != lower) && !n->isSampledAncestorTip() ) )
                {
                    std::string label = "age_" + StringUtilities::toString( (int)(report.node_ages_HPD * 100) ) + "%_HPD";
                    n->addNodeParameter(label, interval);
                }
            }
            else if ( n->isRoot() == false )
            {
                std::string label = "brlen_" + StringUtilities::toString( (int)(report.node_ages_HPD * 100) ) + "%_HPD";
                n->addBranchParameter(label, interval);
            }
        }

    }

    if ( report.node_ages && clock && report.force_positive_branch_lengths )
    {
        enforceNonnegativeBranchLengths( tree.getRoot() );
    }

    if ( report.MAP_parameters )
    {
        mapParameters( tree, verbose );
    }

}


double TreeSummary::cladeProbability( const Clade &c, bool verbose )
{
    summarize( verbose );

    Clade tmp = c;
    tmp.resetTaxonBitset( traces.front()->objectAt(0).getTaxonBitSetMap() );

    return splitFrequency( Split( tmp.getBitRepresentation(), tmp.getMrca(), rooted) );
}


double TreeSummary::computeEntropy( double credible_interval_size, bool stochastic, bool verbose )
{
    std::vector< std::pair<Tree, std::int64_t> > credible_set = getCredibleSetOfTrees(credible_interval_size, stochastic, false);
    int num_taxa = (int)traces.front()->objectAt(0).getTaxa().size();
    double total_samples = sampleSize(true);
    double entropy = 0.0;

    ProgressBar progress = ProgressBar(sampleSize(true));
    
    if (verbose)
    {
        std::stringstream ss;
        ss << "Calculating entropy for a sample of " << credible_set.size() << " unique topologies ..." << std::endl;
        RBOUT( ss.str() );
        progress.start();
    }
    
    for (size_t i = 0; i < credible_set.size(); ++i)
    {
        double freq = credible_set[i].second;
        double p = freq / total_samples;
        entropy += (p * log(p));
        if (verbose) progress.update(i);
    }
    
    if (verbose) progress.finish();

    /* This calculation is directly from AMP / Jeremy's paper. */
    double ln_ntopologies = RbMath::lnFactorial(2 * num_taxa - 5) - RbMath::lnFactorial(num_taxa - 3) - (num_taxa - 3) * RbConstants::LN2;
    entropy += ln_ntopologies;

    return entropy;
}


std::vector<double> TreeSummary::computePairwiseRFDistance( double credible_interval_size, bool stochastic, bool verbose )
{
    std::vector< std::pair<Tree, std::int64_t> > credible_set = getCredibleSetOfTrees(credible_interval_size, stochastic, false);
    std::vector< std::unique_ptr< std::vector<RbBitSet> > > unique_trees_bs(credible_set.size());
    std::vector<size_t> sample_count( credible_set.size() );
    std::vector<double> rf_distances;
    
    ProgressBar progress = ProgressBar(sampleSize(true));
    
    if (verbose)
    {
        std::stringstream ss;
        ss << "Calculating pairwise Robinson-Foulds distances among " << credible_set.size() << " unique topologies ..." << std::endl;
        RBOUT( ss.str() );
        progress.start();
    }
    
    // populate the vector of clade bitsets for all trees
    for (size_t i = 0; i < credible_set.size(); ++i)
    {
        unique_trees_bs[i] = std::make_unique< std::vector<RbBitSet> >();
        credible_set[i].first.getRoot().getAllClades(*(unique_trees_bs[i]), credible_set[i].first.getNumberOfTips(), true);
        sample_count[i] = credible_set[i].second;
    }
    
    for (size_t i = 0; i < credible_set.size(); ++i)
    {
        // The unique tree occurs sample_count[i] times.
        // Here we are treating them as coming in one continuous block, which is annoying.
        
        for (size_t rep = 0; rep < sample_count[i]; rep++)
        {
            // first we need to compare the tree to subsequent copies of itself
            for (size_t k = rep+1; k < sample_count[i]; ++k)
            {
                rf_distances.push_back( 0.0 );
            }
            
            // then we compare it to copies of other trees
            for (size_t j = i+1; j < credible_set.size(); ++j)
            {
                double rf = TreeUtilities::computeRobinsonFouldDistance(*(unique_trees_bs[i]), *(unique_trees_bs[j]), true);
                
                for (size_t k = 0; k < sample_count[j]; ++k)
                {
                    rf_distances.push_back( rf );
                }
            }
        }
        
        if (verbose) progress.update(i);
    }
    
    if (verbose) progress.finish();

    return rf_distances;
}


std::vector<double> TreeSummary::computeTreeLengths( void )
{
    std::vector<double> tree_lengths;

    for(auto& trace: traces)
    {
        for (size_t i = trace->getBurnin(); i < trace->size(); ++i)
        {
            const Tree &tree = trace->objectAt(i);
            tree_lengths.push_back( tree.getTreeLength() );

        }
    }

    return tree_lengths;
}


/* This function, and all other functions that call it, takes a "stochastic" argument that determines what we should do
 * when we cannot obtain a cumulative probability precisely equal to the value passed via the "credible_interval_size"
 * argument. Consider the example worked out by Huelsenbeck & Rannala (2004; Syst. Biol. 53(6): 904--913): the cumulative
 * posterior probability of the first n trees is 0.90, and the probability of the (n+1)-th tree is 0.07. If "stochastic"
 * is true, we will include this tree with include_prob = (0.95 - 0.90) / 0.07 = 0.714. If "stochastic" is false, we will
 * always include it. Note that setting "stochastic" to true can produce an empty credible set if the probability of the
 * most probable tree exceeds "credible_interval_size". E.g., if a single tree accounts for 98% of the total probability
 * mass, we will produce an empty credible set with a probability of 1 - 0.95 / 0.98 = 0.031.
 */
std::vector< std::pair<Tree, std::int64_t> > TreeSummary::getCredibleSetOfTrees( double credible_interval_size, bool stochastic, bool verbose )
{
    summarize( false );

    std::vector< std::pair<Tree, std::int64_t> > unique_trees;
    NewickConverter converter;
    double total_prob = 0;
    double total_samples = sampleSize(true);
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    ProgressBar progress = ProgressBar(sampleSize(true));
    size_t counter = 0;
    
    if (verbose)
    {
        std::stringstream ss;
        ss << "Constructing a " << 100 * credible_interval_size << "% credible set of trees ..." << std::endl;
        RBOUT( ss.str() );
        progress.start();
    }
    
    for (auto& [newick, count]: tree_samples | views::reverse)
    {
        if (verbose) progress.update(counter);
        counter++;
        
        double freq = count;
        double p = freq/total_samples;
        double include_prob = (credible_interval_size - total_prob) / p;
        
        if ( not stochastic or include_prob > rng->uniform01() )
        {
            Tree* current_tree = converter.convertFromNewick( newick );
            current_tree->suppressOutdegreeOneNodes(true);
            std::pair<Tree, std::int64_t> tree_to_add( *current_tree, count );
            unique_trees.push_back( tree_to_add );
            delete current_tree;
        }
        
        total_prob += p;
        
        if ( total_prob >= credible_interval_size )
        {
            break;
        }
    }
    
    if (verbose) progress.finish();

    return unique_trees;
}


std::int64_t TreeSummary::getTopologyCount(const RevBayesCore::Tree &tree, bool verbose)
{
    summarize( verbose );

    Tree t = tree;

    if ( t.isRooted() == false && rooted == false )
    {
        if( outgroup )
        {
            t.reroot( *outgroup, false, true );
        }
        else
        {
            std::vector<std::string> tip_names = traces.front()->objectAt(0).getTipNames();
            std::sort(tip_names.begin(),tip_names.end());
            const std::string& my_outgroup = tip_names[0];
            t.reroot( my_outgroup, false, true );
        }
    }

    std::string newick = t.getPlainNewickRepresentation();

    auto iter = tree_counts.find(newick);
    if (iter == tree_counts.end())
        return 0;
    else
        return iter->second;
}


double TreeSummary::getTopologyFrequency(const RevBayesCore::Tree &tree, bool verbose)
{
    return getTopologyCount(tree, verbose) / (double)sampleSize(true);
}


std::vector<Clade> TreeSummary::getUniqueClades( double min_clade_probability, bool non_trivial_only, bool verbose )
{
    summarize( false );

    std::vector<Clade> unique_clades;
    double total_samples = sampleSize(true);

    std::vector<Taxon> ordered_taxa = traces.front()->objectAt(0).getTaxa();
    VectorUtilities::sort( ordered_taxa );
    size_t num_taxa = ordered_taxa.size();
    
    ProgressBar progress = ProgressBar(sampleSize(true));
    size_t counter = 0;
    
    if (verbose)
    {
        RBOUT("Extracting unique clades ...\n");
        progress.start();
    }

    for (auto& [clade, count]: clade_samples | views::reverse)
    {
        if (verbose) progress.update(counter);
        counter++;
        
        double freq = count;
        double p = freq/total_samples;

        // first we check if this clade is above the minimum level
        if ( p < min_clade_probability )
        {
            break;
        }

        // now let's actually construct the clade
        Clade current_clade(clade.first, ordered_taxa);
        current_clade.setMrca(clade.second);
        
        if (non_trivial_only)
        {
            if ( current_clade.size() == 1 || current_clade.size() == ( rooted ? num_taxa : (num_taxa-1) ) ) continue;
        }

        unique_clades.push_back( current_clade );
    }
    
    if (verbose) progress.finish();

    return unique_clades;
}


std::vector<Tree> TreeSummary::getUniqueTrees( double credible_interval_size, bool stochastic, bool verbose )
{
    std::vector< std::pair<Tree, std::int64_t> > credible_set = getCredibleSetOfTrees(credible_interval_size, stochastic, verbose);
    std::vector<Tree> unique_trees( credible_set.size() );
    
    for (size_t i = 0; i < credible_set.size(); i++)
    {
        unique_trees[i] = credible_set[i].first;
    }

    return unique_trees;
}


bool TreeSummary::isCoveredInInterval(const std::string &v, double credible_interval_size, bool stochastic, bool verbose)
{
    Tree tree;
    tree.initFromString(v);

    if ( tree.isRooted() == false && rooted == false )
    {
        if ( outgroup )
        {
            tree.reroot( *outgroup, false, true );
        }
        else
        {
            std::vector<std::string> tip_names = traces.front()->objectAt(0).getTipNames();
            std::sort(tip_names.begin(),tip_names.end());
            const std::string& this_outgroup = tip_names[0];
            tree.reroot( this_outgroup, false, true );
        }
    }

    return isCoveredInInterval(tree, credible_interval_size, stochastic, verbose);
}


bool TreeSummary::isCoveredInInterval(const Tree &tree, double credible_interval_size, bool stochastic, bool verbose)
{
    summarize( false );
    
    std::vector<Tree> unique_trees = getUniqueTrees(credible_interval_size, stochastic, verbose);
    bool out = false;
    
    for (size_t i = 0; i < unique_trees.size(); i++)
    {
        // comparison of trees compares the strings storing their plain Newick representations
        out |= (unique_trees[i] == tree);
    }

    return out;
}


bool TreeSummary::isClock(void) const
{
    return clock;
}


bool TreeSummary::isDirty(void) const
{
    if (not computed) return true;

    for (auto& trace: traces)
    {
        if (trace->isDirty()) return true;
    }

    return false;
}


double TreeSummary::jointCladeProbability(const RbVector<Clade> &c, bool verbose )
{
    summarize( false );
    
    RbVector<Clade> ref_clades = c;

    size_t num_clades = ref_clades.size();
    for ( size_t i=0; i<num_clades; ++i )
    {
        Clade& tmp = ref_clades[i];
        tmp.resetTaxonBitset( traces.front()->objectAt(0).getTaxonBitSetMap() );
    }
    
    std::vector<std::string> tip_names = traces.front()->objectAt(0).getTipNames();
    std::sort(tip_names.begin(),tip_names.end());
    const std::string& this_outgroup = tip_names[0];

    rooted = traces.front()->objectAt(0).isRooted();

    ProgressBar progress = ProgressBar(sampleSize(true));

    if (verbose)
    {
        std::stringstream ss;
        ss << "Calculating the joint probability of " << c.size() << " clades from " << sampleSize(true) << " trees ..." << std::endl;
        RBOUT( ss.str() );
        progress.start();
    }

    size_t count = 0;
    double num_matches = 0;

    for (auto& trace: traces)
    {
        for (size_t i = trace->getBurnin(); i < trace->size(); ++i)
        {
            if (verbose) progress.update(count);
            count++;

            Tree tree = trace->objectAt(i);

            if ( rooted == false )
            {
                if ( outgroup )
                {
                    tree.reroot( *outgroup, false, true );
                }
                else
                {
                    tree.reroot( this_outgroup, false, true );
                }
            }
            
            bool matches = true;
            for ( size_t i=0; i<num_clades; ++i )
            {
                const Clade& this_clade = ref_clades[i];
                bool contains = tree.getRoot().containsClade( this_clade.getBitRepresentation(), true);
                if ( contains == false )
                {
                    matches = false;
                    break;
                }
                
            }
            
            if ( matches == true )
            {
                ++num_matches;
            }
            
        }
    }

    if (verbose) progress.finish();

    return num_matches / double(count);
}


double TreeSummary::maxdiff( bool verbose )
{
    if (traces.size() <= 1)
    {
        throw RbException("At least 2 traces are required to compute maxdiff");
    }

    std::set<Sample<Split> > splits_union;

    for (auto& trace: traces)
    {
        trace->summarize(verbose);

        splits_union.insert(trace->clade_samples.begin(), trace->clade_samples.end());
    }


    double maxdiff = 0;

    for (auto& [split, count]: splits_union)
    {
        std::vector<double> split_freqs;

        for(auto& trace: traces)
        {
            double freq = trace->splitFrequency(split);

            split_freqs.push_back(freq);
        }

        for(size_t i = 0; i < split_freqs.size(); i++)
        {
            for(size_t j = i+1; j < split_freqs.size(); j++)
            {
                double diff = abs(split_freqs[i] - split_freqs[j]);

                if(diff > maxdiff)
                {
                    maxdiff = diff;
                }
            }
        }
    }

    return maxdiff;
}


Tree* TreeSummary::mapTree( AnnotationReport report, bool verbose )
{
    std::stringstream ss;
    ss << "Compiling maximum a posteriori tree from " << sampleSize(true) << " trees.\n";
    RBOUT(ss.str());

    summarize( verbose );

    // get the tree with the highest posterior probability
    std::string bestNewick = tree_samples.rbegin()->first;
    NewickConverter converter;
    Tree* tmp_best_tree = converter.convertFromNewick( bestNewick );
    tmp_best_tree->suppressOutdegreeOneNodes(true);

    Tree* tmp_tree = NULL;

    if ( clock == true )
    {
        tmp_tree = TreeUtilities::convertTree( *tmp_best_tree );
    }
    else
    {
        tmp_tree = tmp_best_tree->clone();
    }

    delete tmp_best_tree;

    TaxonMap tm = TaxonMap( traces.front()->objectAt(0) );
    tmp_tree->setTaxonIndices( tm );

    report.MAP_parameters = true;
    report.node_ages      = true;
    annotateTree(*tmp_tree, report, false);

    return tmp_tree;
}


Tree* TreeSummary::mccTree( AnnotationReport report, bool verbose )
{
    std::stringstream ss;
    ss << "Compiling maximum clade credibility tree from " << sampleSize(true) << " trees.\n";
    RBOUT(ss.str());

    summarize( verbose );

    Tree* best_tree = NULL;
    std::optional<double> max_cc;

    // find the clade credibility score for each tree
    for (const auto& [newick, count]: tree_samples)
    {
        // find the product of the clade frequencies
        double cc = 0;
        for (auto& [clade, age]: tree_clade_ages.at(newick))
            cc += log( splitFrequency(clade) );

        if (not max_cc or cc > *max_cc)
        {
            max_cc = cc;

            delete best_tree;

            NewickConverter converter;
            Tree* tmp_tree = converter.convertFromNewick( newick );
            tmp_tree->suppressOutdegreeOneNodes(true);
            if ( clock == true )
            {
                best_tree = TreeUtilities::convertTree( *tmp_tree );
            }
            else
            {
                best_tree = tmp_tree->clone();
            }

            TaxonMap tm = TaxonMap( traces.front()->objectAt(0) );
            best_tree->setTaxonIndices( tm );

            delete tmp_tree;
        }
    }

    report.node_ages = true;
    annotateTree(*best_tree, report, false);

    return best_tree;
}


Tree* TreeSummary::mrTree(AnnotationReport report, double cutoff, bool verbose)
{
    if (cutoff < 0.0 || cutoff > 1.0) cutoff = 0.5;

    std::stringstream ss;
    ss << "Compiling majority rule consensus tree (cutoff = " << cutoff << ") from " << sampleSize(true) << " trees.\n";
    RBOUT(ss.str());

    //fill in clades, use all above 50% to resolve the bush with the consensus partitions
    summarize( verbose );        //fills std::vector<Sample<std::string> > cladeSamples, sorts them by descending freq

    //set up variables for consensus tree assembly
    std::vector<std::string> tipNames = traces.front()->objectAt(0).getTipNames();

    //first create a bush
    TopologyNode* root = new TopologyNode(tipNames.size()); //construct root node with index = nb Tips

    for (size_t i = 0; i < tipNames.size(); i++)
    {
        TopologyNode* tipNode = new TopologyNode(tipNames.at(i), i); //Topology node constructor adding tip name and index=taxon nb

        // set the parent-child relationship
        root->addChild(tipNode);
        tipNode->setParent(root);
    }

    //now put the tree together
    Tree* consensusTree = new Tree();
    consensusTree->setRoot(root, true);

    size_t nIndex = tipNames.size();

    double totalSamples = sampleSize(true);

    for (const auto& [clade, count]: clade_samples | views::reverse)
    {
        float cladeFreq = count / totalSamples;
        if (cladeFreq < cutoff)  break;

        //make sure we have an internal node
        size_t clade_size = clade.first.count();
        if (clade_size == 1 || clade_size == tipNames.size())  continue;

        //find parent node
        std::vector<TopologyNode*> children;
        RbBitSet tmp(tipNames.size());
        TopologyNode* parentNode = findParentNode(*root, clade, children, tmp );

        if (parentNode != NULL )
        {
            // skip this node if we've already found a clade compatible with it
            if ( children.size() == parentNode->getNumberOfChildren() ) continue;

            std::vector<TopologyNode*> mrca;

            // find the mrca child if it exists
            if ( clade.second.empty() == false )
            {
                for (size_t i = 0; i < children.size(); i++)
                {
                    if ( children[i]->isTip() && std::find(clade.second.begin(), clade.second.end(), children[i]->getTaxon() ) != clade.second.end() )
                    {
                        mrca.push_back(children[i]);
                    }
                }

                // if we couldn't find the mrca, then this clade is not compatible
                if ( mrca.size() != clade.second.size() )
                {
                    continue;
                }
                else
                {
                    for (size_t i = 0; i < mrca.size(); i++)
                    {
                        mrca[i]->setSampledAncestor(true);
                    }
                }
            }

            nIndex++;   //increment node index
            TopologyNode* intNode = new TopologyNode(nIndex); //Topology node constructor, with proper node index

            // move the children to a new internal node
            for (size_t i = 0; i < children.size(); i++)
            {
                parentNode->removeChild(children[i]);
                intNode->addChild(children[i]);
                children[i]->setParent(intNode);
            }

            intNode->setParent(parentNode);
            parentNode->addChild(intNode);

            // add a mrca child if it exists and there is more than one non-mrca taxa
            if ( mrca.empty() == false && children.size() > 2 )
            {
                TopologyNode* old_parent = parentNode;

                nIndex++;   //increment node index
                parentNode = new TopologyNode(nIndex); //Topology node constructor, with proper node index

                intNode->removeChild(mrca[0]);
                parentNode->addChild(mrca[0]);
                mrca[0]->setParent(parentNode);

                old_parent->removeChild(intNode);
                old_parent->addChild(parentNode);
                parentNode->setParent(old_parent);

                parentNode->addChild(intNode);
                intNode->setParent(parentNode);
            }
        }

        root->setTree(consensusTree);
    }

    //now put the tree together
    consensusTree->setRoot(root, true);

    report.conditional_clade_ages  = false;
    report.conditional_clade_probs = false;
    report.conditional_tree_ages   = false;
    report.node_ages               = true;
    annotateTree(*consensusTree, report, false);

    return consensusTree;
}


void TreeSummary::printCladeSummary(std::ostream &o, double min_clade_probability, bool verbose)
{
    summarize( verbose );

    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(4);

    o << std::endl;
    o << "=========================================" << std::endl;
    o << "Printing Posterior Distribution of Clades" << std::endl;
    o << "=========================================" << std::endl;
    o << std::endl;

    // now the printing
    std::string s = "Samples";
    StringUtilities::fillWithSpaces(s, 16, true);
    o << "\n" << s;
    s = "Posterior";
    StringUtilities::fillWithSpaces(s, 16, true);
    o << s;
    /*s = "ESS";
     StringUtilities::fillWithSpaces(s, 16, true);
     o << s;*/
    s = "Clade";
    StringUtilities::fillWithSpaces(s, 16, true);
    o << s;
    o << std::endl;
    o << "--------------------------------------------------------------" << std::endl;

    double totalSamples = sampleSize(true);

    std::vector<Taxon> ordered_taxa = traces.front()->objectAt(0).getTaxa();
    VectorUtilities::sort( ordered_taxa );

    for (auto& [clade, count]: clade_samples | views::reverse)
    {
        Clade c(clade.first, ordered_taxa);
        c.setMrca(clade.second);

        if ( c.size() == 1 ) continue;

        double freq = count;
        double p = freq/totalSamples;


        if ( p < min_clade_probability )
        {
            break;
        }

        ss.str(std::string());
        ss << freq;
        s = ss.str();
        StringUtilities::fillWithSpaces(s, 16, true);
        o << s;

        ss.str(std::string());
        ss << p;
        s = ss.str();
        StringUtilities::fillWithSpaces(s, 16, true);
        o << s;

        /*ss.str(std::string());
         if ( it->getFrequency() <  totalSamples - burnin && it->getFrequency() > 0 )
         {
         ss << it->getEss();
         }
         else
         {
         ss << " - ";

         }
         s = ss.str();
         StringUtilities::fillWithSpaces(s, 16, true);
         o << s;*/

        o << c;
        o << std::endl;

    }

    o << std::endl;
    o << std::endl;

}


void TreeSummary::printTreeSummary(std::ostream &o, double credible_interval_size, bool stochastic, bool verbose)
{
    summarize( verbose );

    std::stringstream ss;
    ss << std::fixed;
    ss << std::setprecision(4);

    o << std::endl;
    o << "========================================" << std::endl;
    o << "Printing Posterior Distribution of Trees" << std::endl;
    o << "========================================" << std::endl;
    o << std::endl;

    // now the printing
    std::string s = "Cum. Prob.";
    StringUtilities::fillWithSpaces(s, 16, true);
    o << s;
    s = "Samples";
    StringUtilities::fillWithSpaces(s, 16, true);
    o << s;
    s = "Posterior";
    StringUtilities::fillWithSpaces(s, 16, true);
    o << s;
    /*s = "ESS";
     StringUtilities::fillWithSpaces(s, 16, true);
     o << s;*/
    s = "Tree";
    StringUtilities::fillWithSpaces(s, 16, true);
    o << s;
    o << std::endl;
    o << "----------------------------------------------------------------" << std::endl;
    
    double total_prob = 0;
    double total_samples = sampleSize(true);
    RandomNumberGenerator *rng = GLOBAL_RNG;
    
    for (auto& [newick, count]: tree_samples | views::reverse)
    {
        double freq = count;
        double p = freq/total_samples;
        double include_prob = (credible_interval_size - total_prob) / p;
        
        if ( not stochastic or include_prob > rng->uniform01() )
        {
            ss.str(std::string());
            ss << total_prob + p;
            s = ss.str();
            StringUtilities::fillWithSpaces(s, 16, true);
            o << s;

            ss.str(std::string());
            ss << freq;
            s = ss.str();
            StringUtilities::fillWithSpaces(s, 16, true);
            o << s;

            ss.str(std::string());
            ss << p;
            s = ss.str();
            StringUtilities::fillWithSpaces(s, 16, true);
            o << s;

            /* ss.str(std::string());
             ss << it->getEss();
             s = ss.str();
             StringUtilities::fillWithSpaces(s, 16, true);
             o << s; */

            o << newick;
            o << std::endl;
        }
        
        total_prob += p;

        if ( total_prob >= credible_interval_size )
        {
            break;
        }
    }

    o << std::endl;
    o << std::endl;
}


std::int64_t TreeSummary::sampleSize(bool post) const
{
    std::int64_t total = 0;

    for(auto& trace: traces)
    {
        total += trace->size(post);
    }

    return total;
}


void TreeSummary::setOutgroup(const RevBayesCore::Clade &c)
{
    outgroup = c;
}


TreeSummary::Split TreeSummary::collectTreeSample(const TopologyNode& n, RbBitSet& intaxa, std::string newick, std::map<Split, std::int64_t>& cladeCountMap)
{
    double age = (clock ? n.getAge() : n.getBranchLength() );

    std::vector<Split> child_splits;

    RbBitSet taxa(intaxa.size());
    std::set<Taxon> mrca;

    if ( n.isTip() )
    {
        n.getTaxa(taxa);

        if ( rooted && n.isSampledAncestorTip() )
        {
            sampled_ancestor_counts[n.getTaxon()]++;

            mrca.insert( n.getTaxon() );
        }
    }
    else
    {
        for (size_t i = 0; i < n.getNumberOfChildren(); i++)
        {
            const TopologyNode &child_node = n.getChild(i);

            child_splits.push_back( collectTreeSample(child_node, taxa, newick, cladeCountMap) );

            if ( rooted && child_node.isSampledAncestorTip() )
            {
                mrca.insert(child_node.getTaxon());
            }
        }
    }

    intaxa |= taxa;

    Split parent_split(taxa, mrca, rooted);

    if ( taxa.size() > 0 )
    {
        // store the age for this split
        clade_ages[parent_split].push_back( age );

        // increment split count
        cladeCountMap[parent_split]++;

        // add conditional clade ages
        for (auto& child_split: child_splits)
        {
            // inserts new entries if doesn't already exist
            conditional_clade_ages[parent_split][child_split].push_back( clade_ages[child_split].back() );
        }

        // store the age for this split, conditional on the tree topology
        tree_clade_ages[newick][parent_split].push_back( age );
    }

    return parent_split;
}


void TreeSummary::enforceNonnegativeBranchLengths(TopologyNode& node) const
{
    std::vector<TopologyNode*> children = node.getChildren();

    double minimum_branch_length = 1e-6;
    for (size_t i = 0; i < children.size(); i++)
    {
        if (children[i]->getAge() > node.getAge())
        {
            children[i]->setAge( node.getAge() - minimum_branch_length );
        }
        enforceNonnegativeBranchLengths( *children[i] );
    }
}


TopologyNode* TreeSummary::findParentNode(TopologyNode& n, const Split& split, std::vector<TopologyNode*>& children, RbBitSet& child_b ) const
{
    size_t num_taxa = child_b.size();

    RbBitSet node( num_taxa );
    n.getTaxa(node);

    RbBitSet clade = split.first;

    RbBitSet mask  = node | clade;

    bool compatible = (mask == node);
    bool ischild      = (mask == clade);

    Split c = split;
    // check if the flipped unrooted split is compatible
    if ( !rooted && !compatible && !ischild)
    {
        RbBitSet clade_flip = ~clade;
        mask  = node | clade_flip;

        compatible = (mask == node);

        if ( compatible )
        {
            c.first = clade_flip;
        }
    }

    TopologyNode* parent = NULL;

    if (compatible)
    {
        parent = &n;

        std::vector<TopologyNode*> new_children;

        // keep track of which taxa we found in the children
        RbBitSet child_mask( num_taxa );

        for (auto& old_child: n.getChildren())
        {
            RbBitSet child_bset(clade.size());

            TopologyNode* child = findParentNode(*old_child, c, new_children, child_bset );

            // add this child to the mask
            child_mask = (child_bset | child_mask);

            // check if child is a compatible parent
            if (child != NULL)
            {
                parent = child;
                break;
            }
        }

        children = std::move(new_children);

        // check that we found all the children
        if ( parent == &n && child_mask != c.first && !n.isTip())
        {
            parent = NULL;
        }
    }
    else if (ischild)
    {
        child_b = node;
        children.push_back(&n);
    }

    return parent;
}


/*
 * this method calculates the MAP ancestral character states for the nodes on the input_tree
 */
void TreeSummary::mapContinuous(Tree &tree, const std::string &n, size_t paramIndex, double hpd, bool isNodeParameter, bool verbose  ) const
{

    // 2-d vectors to keep the data (posteriors and states) of the inputTree nodes: [node][data]
    const std::vector<TopologyNode*> &summary_nodes = tree.getNodes();
    std::vector<std::vector<double> > samples(summary_nodes.size(),std::vector<double>());

    // flag if only interior nodes are used
    bool interior_only  = false;
    bool tips_checked   = false;
    bool root_checked   = false;
    bool use_root       = true;

    // get the reference newick representation
    std::string reference_newick = tree.getPlainNewickRepresentation();

    // start the progress bar
    ProgressBar progress = ProgressBar(sampleSize(true));
    if ( verbose == true )
    {
        progress.start();
    }

    // for fast access, get the newick strings for all clades
    std::vector<std::string> summary_newick;
    size_t num_nodes = summary_nodes.size();
    for (size_t j = 0; j < num_nodes; ++j)
    {
        TopologyNode *node = summary_nodes[j];
        summary_newick.push_back( node->computePlainNewick() );
    }

    size_t count = 0;

    for(auto& trace: traces)
    {
        // loop through all trees in tree trace
        for (size_t i = trace->getBurnin(); i < trace->size(); ++i)
        {
            if ( verbose == true )
            {
                progress.update( count );
                count++;
            }

            const Tree &sample_tree = trace->objectAt( i );
            const TopologyNode& sample_root = sample_tree.getRoot();

            // create a map from newick strings to clade indices
            // we also get the newick representation of this sample
            std::map<std::string,size_t> sample_clade_indices;
            std::string sample_newick = sample_root.fillCladeIndices(sample_clade_indices);

            // compare if this tree is the same as the reference tree
            bool same_tree = ( reference_newick == sample_newick );

            // loop through all nodes in inputTree
            size_t num_nodes = summary_nodes.size();
            for (size_t j = 0; j < num_nodes; ++j)
            {
                TopologyNode *node = summary_nodes[j];
                if ( node->isTip() == true )
                {
                    if ( tips_checked == false )
                    {

                        tips_checked = true;
                        size_t sample_clade_index = sample_clade_indices[ summary_newick[j] ];

                        const TopologyNode &sample_node = sample_tree.getNode( sample_clade_index );

                        std::vector<std::string> params;
                        if ( isNodeParameter == true )
                        {
                            params = sample_node.getNodeParameters();
                        }
                        else
                        {
                            params = sample_node.getBranchParameters();
                        }

                        // check if this parameter exists
                        if ( params.size() > paramIndex )
                        {

                            std::string tmp = params[paramIndex];
                            if ( tmp[0] == '&')
                            {
                                tmp = tmp.substr(1,tmp.size());
                            }
                            std::vector<std::string> pair;
                            StringUtilities::stringSplit(tmp, "=", pair);

                            // check if this parameter has the correct name
                            interior_only = (pair[0] != n);
                        }
                        else
                        {
                            interior_only = true;
                        }


                    }

                    if ( interior_only == true )
                    {
                        continue;
                    }
                }

                if ( node->isRoot() == true )
                {
                    if ( root_checked == false )
                    {

                        root_checked = true;

                        const TopologyNode &sample_node = sample_tree.getRoot();

                        std::vector<std::string> params;
                        if ( isNodeParameter == true )
                        {
                            params = sample_node.getNodeParameters();
                        }
                        else
                        {
                            params = sample_node.getBranchParameters();
                        }

                        // check if this parameter exists
                        if ( params.size() > paramIndex )
                        {

                            std::string tmp = params[paramIndex];
                            if ( tmp[0] == '&')
                            {
                                tmp = tmp.substr(1,tmp.size());
                            }
                            std::vector<std::string> pair;
                            StringUtilities::stringSplit(tmp, "=", pair);

                            // check if this parameter has the correct name
                            use_root = pair[0] == n;
                        }
                        else
                        {
                            use_root = false;
                        }


                    }

                    if ( use_root == false )
                    {
                        continue;
                    }

                }


                if ( same_tree == true || sample_clade_indices.find( summary_newick[j] ) != sample_clade_indices.end() )
                {
                    // if the inputTree node is also in the sample tree
                    // we get the ancestral character state from the ancestral state trace
                    size_t sample_clade_index = sample_clade_indices[ summary_newick[j] ];
                    const TopologyNode &sample_node = sample_tree.getNode( sample_clade_index );

                    std::vector<std::string> params;
                    if ( isNodeParameter == true )
                    {
                        params = sample_node.getNodeParameters();
                    }
                    else
                    {
                        params = sample_node.getBranchParameters();
                    }

                    // check if this parameter exists
                    if ( params.size() <= paramIndex )
                    {
                        throw RbException() << "Too few parameters for this tree during the tree annotation. Problematic tree: " << sample_tree.getNewickRepresentation();
                    }

                    std::string tmp = params[paramIndex];
                    if ( tmp[0] == '&')
                    {
                        tmp = tmp.substr(1,tmp.size());
                    }
                    std::vector<std::string> pair;
                    StringUtilities::stringSplit(tmp, "=", pair);

                    // check if this parameter has the correct name
                    if ( pair[0] != n )
                    {

                        throw RbException("The parameter for this tree doesn't match during the tree annotation.");
                    }

                    double state = atof(pair[1].c_str());

                    std::vector<double> &entries = samples[j];
                    entries.push_back( state );

                } // end if the sampled tree contained this clade

            } // end loop over all nodes in the tree

        } // end loop over each iteration in the trace

    } // end loop over each trace


    // end the progress bar
    if ( verbose == true )
    {
        progress.finish();
    }

    std::vector<double> posteriors;
    for (int idx = 0; idx < num_nodes; ++idx)
    {

        TopologyNode &node = *summary_nodes[idx];
        if ( ( node.isTip() == false || interior_only == false ) && ( node.isRoot() == false || use_root == true ) )
        {

            // collect the samples
            std::vector<double> stateSamples = samples[idx];

            // sort the samples by frequency
            sort(stateSamples.begin(), stateSamples.end());


            size_t interval_start = ((1.0-hpd)/2.0) * stateSamples.size();
            size_t interval_median = 0.5 * stateSamples.size();
            size_t interval_end = (1.0-(1.0-hpd)/2.0) * stateSamples.size();
            interval_end = (interval_end >= stateSamples.size() ? stateSamples.size()-1 : interval_end);
            double lower = stateSamples[interval_start];
            double median = stateSamples[interval_median];
            double upper = stateSamples[interval_end];

            // make node age annotation
            std::string param = "{" + StringUtilities::toString(lower)
            + "," + StringUtilities::toString(upper) + "}";

            if ( isNodeParameter == true )
            {
                // make parameter string for this node
                node.addNodeParameter(n+"_range",param);

                // make parameter string for this node
                node.addNodeParameter(n,median);
            }
            else
            {

                // make parameter string for this node
                node.addBranchParameter(n+"_range",param);

                // make parameter string for this node
                node.addBranchParameter(n,median);

            }

        }

    }

}


/*
 * this method calculates the MAP ancestral character states for the nodes on the input_tree
 */
void TreeSummary::mapDiscrete(Tree &tree, const std::string &n, size_t paramIndex, size_t num, bool isNodeParameter, bool verbose ) const
{

    // 2-d vectors to keep the data (posteriors and states) of the inputTree nodes: [node][data]
    const std::vector<TopologyNode*> &summary_nodes = tree.getNodes();
    //std::vector<std::map<std::string, Sample<std::string> > > stateAbsencePresence(summary_nodes.size(), std::map<std::string, Sample<std::string> >());
    std::vector<std::map<std::string, std::int64_t> > state_counts(summary_nodes.size(), std::map<std::string, std::int64_t>());

    bool interiorOnly = true;
    bool tipsChecked = false;
    //    bool useRoot = true;

    size_t total_size = 0;

    for(auto& trace: traces)
    {
        total_size += trace->size(true);

        // loop through all trees in tree trace
        for (size_t iteration = trace->getBurnin(); iteration < trace->size(); ++iteration)
        {
            const Tree &sample_tree = trace->objectAt( iteration );
            const TopologyNode& sample_root = sample_tree.getRoot();

            // loop through all nodes in inputTree
            for (size_t node_index = 0; node_index < summary_nodes.size(); ++node_index)
            {
                TopologyNode *node = summary_nodes[node_index];

                if ( node->isTip() == true )
                {
                    if ( tipsChecked == false )
                    {
                        tipsChecked = true;
                        size_t sample_clade_index = sample_root.getCladeIndex( node );

                        const TopologyNode &sample_node = sample_tree.getNode( sample_clade_index );

                        std::vector<std::string> params;
                        if ( isNodeParameter == true )
                        {
                            params = sample_node.getNodeParameters();
                        }
                        else
                        {
                            params = sample_node.getBranchParameters();
                        }

                        // check if this parameter exists
                        if ( params.size() > paramIndex )
                        {

                            std::string tmp = params[paramIndex];
                            if ( tmp[0] == '&')
                            {
                                tmp = tmp.substr(1,tmp.size());
                            }
                            std::vector<std::string> pair;
                            StringUtilities::stringSplit(tmp, "=", pair);

                            // check if this parameter has the correct name
                            interiorOnly = pair[0] != n;
                        }


                    }

                    if ( interiorOnly == true )
                    {
                        continue;
                    }
                }

                if ( sample_root.containsClade(node, true) )
                {
                    // if the inputTree node is also in the sample tree
                    // we get the ancestral character state from the ancestral state trace
                    size_t sample_clade_index = sample_root.getCladeIndex( node );

                    const TopologyNode &sample_node = sample_tree.getNode( sample_clade_index );

                    std::vector<std::string> params;
                    if ( isNodeParameter == true )
                    {
                        params = sample_node.getNodeParameters();
                    }
                    else
                    {
                        params = sample_node.getBranchParameters();
                    }

                    // check if this parameter exists
                    if ( params.size() <= paramIndex )
                    {
                        if ( sample_node.isRoot() == true )
                        {
                            continue;
                        }
                        else
                        {
                            throw RbException() << "Too few parameters for this tree during the tree annotation. Problematic tree: " << sample_tree.getNewickRepresentation();
                        }

                    }

                    std::string tmp = params[paramIndex];
                    if ( tmp[0] == '&')
                    {
                        tmp = tmp.substr(1,tmp.size());
                    }
                    std::vector<std::string> pair;
                    StringUtilities::stringSplit(tmp, "=", pair);

                    // check if this parameter has the correct name
                    if ( pair[0] != n )
                    {
                        throw RbException("The parameter for this tree doesn't match during the tree annotation.");
                    }

                    const std::string &state = pair[1];

                    state_counts[node_index][state]++;

                } // end if the sampled tree contained this clade

            } // end loop over all nodes in the tree

        } // end loop over each iteration in the trace

    } // end loop over each trace


    std::vector<double> posteriors;
    for (int i = 0; i < summary_nodes.size(); i++)
    {

        TopologyNode &node = *summary_nodes[i];
        if ( node.isTip() && interiorOnly == true )
        {
            // make parameter string for this node
            if ( isNodeParameter == true )
            {
                node.addNodeParameter(n,"{}");
            }
            else
            {
                node.addBranchParameter(n,"{}");
            }
        }
        else
        {

            // collect the samples
            std::set<Sample<std::string> > stateSamples;
            for (auto& [state,count]: state_counts[i])
            {
                stateSamples.insert( {state, count} );
            }

            double total_node_pp = 0.0;
            std::string final_state = "{";
            int i=0;
            for (auto& [state, count]: stateSamples)
            {
                if ( total_node_pp > 0.9999 ) continue;

                if (i > 0)
                    final_state += ",";
                i++;

                double pp = count / total_size;
                final_state += state + "=" + StringUtilities::toString(pp);
                total_node_pp += pp;

            }

            final_state += "}";

            // make parameter string for this node
            if ( isNodeParameter == true )
            {
                node.addNodeParameter(n,final_state);
            }
            else
            {
                node.addBranchParameter(n,final_state);
            }
        }
    }

}


// annotate the MAP node/branch parameters
void TreeSummary::mapParameters( Tree &tree, bool verbose ) const
{

    const Tree& sample_tree = traces.front()->objectAt( 0 );

    // first we annotate the node parameters
    // we need an internal node because the root might not have all parameter (e.g. rates)
    // and the tips might neither have all parameters
    const TopologyNode *n = &sample_tree.getRoot().getChild( 0 );
    if ( n->isTip() == true )
    {
        n = &sample_tree.getRoot().getChild( 1 );
    }
    const std::vector<std::string> &nodeParameters = n->getNodeParameters();
    for (size_t i = 0; i < nodeParameters.size(); ++i)
    {

        std::string tmp = nodeParameters[i];
        if ( tmp[0] == '&')
        {
            tmp = tmp.substr(1,tmp.size());
        }
        std::vector<std::string> pair;
        StringUtilities::stringSplit(tmp, "=", pair);

        if ( pair[0] == "index" ) continue;

        if ( StringUtilities::isNumber( pair[1] ) && !StringUtilities::isIntegerNumber( pair[1] ) )
        {
            mapContinuous(tree, pair[0], i, 0.95, true, verbose);
        }
        else
        {
            mapDiscrete(tree, pair[0], i, 3, true, verbose);
        }

    }

    // then we annotate the branch parameters
    const std::vector<std::string> &leftBranchParameters = sample_tree.getRoot().getChild(0).getBranchParameters();
    const std::vector<std::string> &rightBranchParameters = sample_tree.getRoot().getChild(1).getBranchParameters();

    std::vector<std::string> branchParameters;
    if ( leftBranchParameters.size() > rightBranchParameters.size() )
    {
        branchParameters = leftBranchParameters;
    }
    else
    {
        branchParameters = rightBranchParameters;
    }

    for (size_t i = 0; i < branchParameters.size(); ++i)
    {

        std::string tmp = branchParameters[i];
        if ( tmp[0] == '&')
        {
            tmp = tmp.substr(1,tmp.size());
        }
        std::vector<std::string> pair;
        StringUtilities::stringSplit(tmp, "=", pair);

        if ( pair[0] == "index" ) continue;

        if ( StringUtilities::isNumber( pair[1] ) )
        {
            mapContinuous(tree, pair[0], i, 0.95, false, verbose);
        }
        else
        {
            mapDiscrete(tree, pair[0], i, 3, false, verbose);
        }

    }

}


std::int64_t TreeSummary::splitCount(const Split &n) const
{
    auto iter = clade_counts.find(n);

    if (iter == clade_counts.end())
        return 0;
    else
        return iter->second;
}


double TreeSummary::splitFrequency(const Split &n) const
{
    return double(splitCount(n))/sampleSize(true);
}


void TreeSummary::summarize( bool verbose )
{
    if ( not isDirty() )
    {
        if (verbose) RBOUT("Skipping the clade summarization step, as a summary was already produced by a previous method call ...\n");
        return;
    }

    std::vector<std::string> tip_names = traces.front()->objectAt(0).getTipNames();
    std::sort(tip_names.begin(),tip_names.end());
    const std::string& this_outgroup = tip_names[0];

    rooted = traces.front()->objectAt(0).isRooted();

    clade_samples.clear();
    tree_samples.clear();

    sampled_ancestor_counts.clear();

    clade_ages.clear();
    conditional_clade_ages.clear();
    tree_clade_ages.clear();

    clade_counts.clear();
    tree_counts.clear();

    ProgressBar progress = ProgressBar(sampleSize(true));

    if (verbose)
    {
        RBOUT("Summarizing clades ...\n");
        progress.start();
    }

    size_t count = 0;

    for (auto& trace: traces)
    {
        for (size_t i = trace->getBurnin(); i < trace->size(); ++i)
        {
            if (verbose) progress.update(count);
            count++;

            Tree tree = trace->objectAt(i);

            if ( rooted == false )
            {
                if ( outgroup )
                {
                    tree.reroot( *outgroup, false, true );
                }
                else
                {
                    tree.reroot( this_outgroup, false, true );
                }
            }

            std::string newick = tree.getPlainNewickRepresentation();

            tree_counts[newick]++;

            // get the clades for this tree
            RbBitSet b( tree.getNumberOfTips(), false );
            collectTreeSample(tree.getRoot(), b, newick, clade_counts);
        }
    }

    // FIXME: Lots of the loops we do use reverse iterators.
    //
    //        so maybe sort in order of descending frequency instead?

    // FIXME: Storing a second copy of everything just to know the order  wastes memory.
    //
    //        Possibly store the clades (and trees) in a vector, and then make
    //        clade_counts/clade_samples (and tree_counts/tree_samples) just refer to the
    //        index in the vector.
    
    // sort the clade samples in ascending frequency
    for (auto& [clade, count]: clade_counts)
    {
//        if ( it->first.first.count() > 0 )
//        if ( it->first.first.count() > 0 && it->first.first.count() < (num_taxa-1) )
        {
            clade_samples.insert( Sample<Split>(clade, count) );
        }

    }

    // sort the tree samples in ascending frequency
    for (auto& [tree, count]: tree_counts)
    {
        tree_samples.insert( {tree, count} );
    }

    // finish progress bar
    if (verbose) progress.finish();

    // Mark all the traces clean.
    for (auto& trace: traces)
    {
        trace->setDirty(false);
    }

    // And record that the summary statistics are computed.
    computed = true;

    // FIXME - Why are we marking input traces dirty?
    //         This isn't really a property of the input trace,
    //         but a property of the TreeSummary and its summary statistics.
}
