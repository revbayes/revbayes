#include "Transforms.h"
#include "TypedDagNode.h"

using std::optional;
using std::vector;

namespace Transforms
{
    using namespace RevBayesCore;
    // f(x) = x + l

    optional<double> add_transform(const vector<const DagNode*>& params, double x)
    {
        auto l = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return x + l->getValue();
    }

    optional<double> add_inverse(const vector<const DagNode*>& params, double x)
    {
        auto l = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return x - l->getValue();
    }

    optional<double> log_add_prime(const vector<const DagNode*>& params, double x)
    {
        return 0;
    }

    // f(x) = x * l

    optional<double> mul_transform(const vector<const DagNode*>& params, double x)
    {
        auto l = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return x * l->getValue();
    }

    optional<double> mul_inverse(const vector<const DagNode*>& params, double x)
    {
        auto l = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return x / l->getValue();
    }

    optional<double> log_mul_prime(const vector<const DagNode*>& params, double x)
    {
        auto l = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return log(std::abs(l->getValue()));
    };

    // f(x) = exp(x)
    optional<double> exp_transform(const vector<const DagNode*>& params, double x)
    {
        return exp(x);
    }

    optional<double> exp_inverse(const vector<const DagNode*>& params, double x)
    {
        if (x > 0)
            return log(x);
        else
            return {}; // out of range
    }

    optional<double> log_exp_prime(const vector<const DagNode*>& params, double x)
    {
        // y = exp(x)
        // dy/dx = exp(x)
        // log(dy/dx) = x;

        return x;
    }
    
    // f(x) = log(x)
    optional<double> log_transform(const vector<const DagNode*>& params, double x)
    {
        if (x > 0)
            return log(x);
        else
            return {}; // out of range
    }

    optional<double> log_inverse(const vector<const DagNode*>& params, double x)
    {
        return exp(x);
    }

    optional<double> log_log_prime(const vector<const DagNode*>& params, double x)
    {
        // y = log(x)
        // dy/dx = 1/x
        // log(dy/dx) = -log(x)

        if (x > 0)
            return -log(x);
        else
            return {}; // out of range
    }


    // f(x) = logit(x) = log(x/(1-x))

    // This can handle x=exp(-1000), but not x=1 - exp(-1000)
    optional<double> logit_transform(const vector<const DagNode*>& params, double x)
    {
        if (x > 0 and x < 1)
        {
            // log( x / (1-x) )
            return log(x) - log1p(-x);
        }
        else
            return {}; // out of range
    }

    optional<double> logit_inverse(const vector<const DagNode*>& params, double x)
    {
        // Two forms of the same expression to avoid overflow with exp(larg number)
        if (x < 0)
            return exp(x)/(1+exp(x));
        else
            return 1/(exp(-x)+1);
    }

    optional<double> log_logit_prime(const vector<const DagNode*>& params, double x)
    {
        // y = log(x/(1-x))
        // dy/dx = 1/(x/(1-x)) * ((1-x)(1) - x(-1))/((1-x)**2)
        //       = (1-x)/x     * 1/((1-x)**2)
        //       = 1/(x * (1-x) )
        // logit(dy/dx) = -log(x) - log(1-x)

        if (x > 0 and x < 1)
            return -log(x) - log1p(-x);
        else
            return {}; // out of range
    }

    // f(x) = x - sec

    optional<double> sub1_transform(const vector<const DagNode*>& params, double x)
    {
        auto sec = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return x - sec->getValue();
    }

    optional<double> sub1_inverse(const vector<const DagNode*>& params, double x)
    {
        auto sec = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return x + sec->getValue();
    }

    optional<double> log_sub1_prime(const vector<const DagNode*>& params, double x)
    {
        return 0;
    }

    // f(x) = fir - x

    optional<double> sub2_transform(const vector<const DagNode*>& params, double x)
    {
        auto fir = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return fir->getValue() - x; 
    }

    optional<double> sub2_inverse(const vector<const DagNode*>& params, double x)
    {
        auto fir = dynamic_cast<const TypedDagNode<double>*>(params[0]);
        return fir->getValue() - x; 
    }

    optional<double> log_sub2_prime(const vector<const DagNode*>& params, double x)
    {
        return 0;
    }

    // f(x) = exp(x)/(1+exp(x)) -- the inverse of the logit transform.
    optional<double> invlogit_transform(const vector<const DagNode*>& params, double x)
    {
        // Two forms of the same expression to avoid overflow with exp(larg number)
        if (x < 0)
            return exp(x)/(1+exp(x));
        else
            return 1/(exp(-x)+1);
    }

    optional<double> invlogit_inverse(const vector<const DagNode*>& params, double x)
    {
        if (x > 0 and x < 1)
        {
            // log( x / (1-x) )
            return log(x) - log1p(-x);
        }
        else
            return {}; // out of range
    }

    optional<double> log_invlogit_prime(const vector<const DagNode*>& params, double x)
    {
        if (x > 0)
        {
            // y = 1/(exp(-x) + 1)
            // dy/dx = -1/((exp(-x)+1)^2) * exp(-x) * -1
            //       = exp(-x) / (1 + exp(-x))^2
            // log(dy/dx) = -x - 2*log(1+exp(-x))

            return - x - 2*log1p(exp(-x));
        }
        else
        {
            // y = exp(x)/(exp(x) + 1)
            // dy/dx = exp(x) / (1+exp(x))^2
            // log(dy/dx) = x  - 2log1p(exp(x))
    
            return   x - 2*log1p(exp(x));
        }
    }


}
