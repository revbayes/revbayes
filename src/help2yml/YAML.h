#ifndef YAML_H
#define YAML_H

#include <string>
#include <vector>

namespace YAML
{
    class Element
    {
        public:

        Element( std::string v ) : value(v) {};
        virtual ~Element() {};

        std::string value;
    };

    class Node;

    class Attribute : public Element
    {
        public:

        Attribute( std::string v );
        Attribute( std::string v, std::vector<Element*> it );

        Attribute( std::string v, std::string it );
        Attribute( std::string v, std::vector<std::string> it );

        ~Attribute();

        std::vector<Node*> children;
    };

    class Node
    {
        public:

        Node();
        Node( const std::vector<Element*>& );
        ~Node();

        std::string render(size_t depth = 0);
        std::string format(std::string value, size_t depth = 0);

        std::vector<Element*> elements;

    };
}

#endif

