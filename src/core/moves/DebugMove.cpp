#include <iomanip>
#include <iostream>
#include <range/v3/all.hpp> // for ranges::views

#include "DebugMove.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbConstants.h"

using namespace RevBayesCore;

namespace views = ranges::views;

std::map<const DagNode*, double> getNodePrs(const std::vector<DagNode*>& nodes, const RbOrderedSet<DagNode*>& affected_nodes)
{
    std::map<const DagNode*, double> Prs;
    for(auto node: views::concat(nodes, affected_nodes))
    {
        double lnPr = RbConstants::Double::nan;
        try
        {
            lnPr = node->getLnProbability();
        }
        catch (RbException& e)
        {
            if (e.getExceptionType() != RbException::MATH_ERROR)
                throw;
        }
        Prs.insert({node, lnPr});
    }
    return Prs;
}

void compareNodePrs(const std::string& name, const std::map<const DagNode*, double>& pdfs1, const std::map<const DagNode*, double>& pdfs2, const std::string& msg)
{
    RbException E;
    E<<std::setprecision(err_precision)<<"Executing "<<name<<": "<<msg<<"!\n";
    bool err = false;
    bool weird = false;
    for(auto& [node,pr1]: pdfs1)
    {
	auto pr2 = pdfs2.at(node);

	// If they are equal then there is no difference.
	if (pr1 == pr2 or (std::isnan(pr1) and std::isnan(pr2)))
        {
            // But if both are NaN or -Inf then its kind of weird!
            // This could indicate that a previous move created a -Inf situation.
            // In the future it could indicate that we haven't yet moved from Pr=0 to Pr>0.
            if (not std::isfinite(pr1))
            {
                std::cerr<<"    WEIRD: "<<node->getName()<<": "<<pr1<<"\n";
                weird = true;
            }
            continue;
        }

	// Be a bit careful about computing a relative error.
	double abs_err = std::abs(pr1 - pr2);
	double scale = std::min(std::abs(pr1),std::abs(pr2));
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

