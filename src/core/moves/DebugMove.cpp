#include <iomanip>
#include <iostream>
#include <range/v3/all.hpp> // for ranges::views
#include <boost/lexical_cast.hpp> // for lexical_cast

#include "DebugMove.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbConstants.h"
#include "variant.h"

using namespace RevBayesCore;

namespace views = ranges::views;

const LogDensity* ProbOrError::is_prob() const
{
    return to<LogDensity>(*this);
}

const std::string* ProbOrError::is_error() const
{
    return to<std::string>(*this);
}

LogDensity ProbOrError::is_prob_or(LogDensity d) const
{
    if (auto p = is_prob())
        return *p;
    else
        return d;
}

LogDensity ProbOrError::is_prob_or_nan() const
{
    return is_prob_or( std::nan("1") );
}

std::ostream& operator<<(std::ostream& o, const ProbOrError& pe)
{
    if (auto d = pe.is_prob())
        o<<*d;
    else if (auto s = pe.is_error())
        o<<"Exception: "<<*s;
    else
        std::abort();

    return o;
}

NodePrMap getNodePrs(const std::vector<DagNode*>& nodes, const RbOrderedSet<DagNode*>& affected_nodes)
{
    NodePrMap Prs;
    for(auto node: views::concat(nodes, affected_nodes))
    {
        ProbOrError lnPr(RbConstants::Double::nan);

        try
        {
            Prs.insert({node, node->getLnProbability()});
        }
        catch (RbException& e)
        {
            if (e.getExceptionType() == RbException::MATH_ERROR)
                Prs.insert({node, e.getMessage()});
            else
                throw;
        }
    }
    return Prs;
}

void compareNodePrs(const std::string& name, const NodePrMap& pdfs1, const NodePrMap& pdfs2, const std::string& msg)
{
    RbException E;
    E<<std::setprecision(err_precision)<<"Executing "<<name<<": "<<msg<<"!\n";
    bool err = false;
    bool weird = false;
    for(auto& [node,PE1]: pdfs1)
    {
	auto PE2 = pdfs2.at(node);

        LogDensity pr1 = PE1.is_prob_or(RbConstants::Double::nan);
        LogDensity pr2 = PE2.is_prob_or(RbConstants::Double::nan);

	// If they are equal then there is no difference.
	if (pr1 == pr2 or (pr1.isnan() and pr2.isnan()))
        {
            // But if both are NaN or -Inf then its kind of weird!
            // This could indicate that a previous move created a -Inf situation.
            // In the future it could indicate that we haven't yet moved from Pr=0 to Pr>0.
            if (not pr1.isfinite())
            {
                std::cerr<<"    WEIRD: "<<node->getName()<<": "<<PE1<<"\n";
                weird = true;
            }
            continue;
        }

	// Be a bit careful about computing a relative error.
	double abs_err = std::abs(double(pr1 - pr2));
	double scale = std::min(std::abs(double(pr1)),std::abs(double(pr2)));
	if (scale < 1 or not RbMath::isAComputableNumber(scale))
	    scale = 1;
	double rel_err = abs_err/scale;

	// Complain if rel_err is NaN.
	if (not (rel_err < rel_err_threshhold))
	{
	    E<<"    "<<node->getName()<<": "<<pr1<<" != "<<pr2<<"    diff = "<<pr1-pr2<<"\n";
	    err = true;
	}
    }

    if (weird) std::cerr<<std::endl;

    if (err) throw E;
}

void showNodePrs(const std::string& stage, const NodePrMap& pdfs)
{
    for(auto& [node, pr]: pdfs)
    {
        std::cerr<<"    "<<stage<<": "<<node->getName()<<":  "<<pr<<"\n";
    }
    std::cerr<<"\n";
}

void showNodePrs(const std::string& stage, const NodePrMap& pdfs, const std::vector<const RevBayesCore::DagNode*>& nodes)
{
    // This version tries to print the nodes that we modify first, and the affected nodes second.
    std::set<const RevBayesCore::DagNode*> nodes_set;

    // First print the PDFS for the nodes that were directly modified.
    for(auto node: nodes)
    {
        auto pr = pdfs.at(node);
        std::cerr<<"    "<<stage<<": "<<node->getName()<<":  "<<pr<<"\n";

        nodes_set.insert(node);
    }

    // These are the pdfs that are NOT in nodes -- presumably the affectedNodes
    for(auto& [node, pr]: pdfs)
    {
        if (not nodes_set.contains(node))
            std::cerr<<"    "<<stage<<": "<<node->getName()<<":  "<<pr<<"\n";
    }

    std::cerr<<"\n";
}

void checkNotTouched(const std::vector<RevBayesCore::DagNode*>& nodes, const RevBayesCore::RbOrderedSet<RevBayesCore::DagNode*>& affected_nodes)
{
    bool ok = true;
    for(auto node: views::concat(nodes, affected_nodes))
    {
        if (node->isTouched())
        {
            std::cerr<<"  node "<<node->getName()<<" is already touched!\n";
            ok = false;
        }
    }
    if (not ok) std::abort();
}
