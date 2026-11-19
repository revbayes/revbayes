#ifndef RevObjectConversionFunction_H
#define RevObjectConversionFunction_H

#include "AbstractModelObject.h"
#include "Cloner.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

#include <memory>

namespace RevLanguage {

template <class T> class ModelObject;

// Preserve the semantics of RevObject::convertTo while making its result depend on a source node.
// A private source copy supplies the Rev type; the parameter supplies the current value after rewiring.
template <class T>
class RevObjectConversionFunction : public RevBayesCore::TypedFunction<T> {
public:
    // Keep the already checked initial conversion and track its source independently of the variable wrapper.
    RevObjectConversionFunction(const AbstractModelObject& source, const ModelObject<T>& converted) :
        RevBayesCore::TypedFunction<T>(copyValue(converted.getValue())),
        sourceNode(source.getDagNode()),
        targetType(converted.getTypeSpec())
    {
        // The prototype retains only the Rev type and a value, not the source's parents.
        std::unique_ptr<AbstractModelObject> prototype(source.clone());
        prototype->makeConstantValue();
        sourceObject = std::move(prototype);
        this->addParameter(sourceNode);
    }

    RevObjectConversionFunction* clone() const override
    {
        return new RevObjectConversionFunction(*this);
    }

    // Convert a private snapshot, so conversion code cannot change or retain the live source node.
    // In particular, element-wise and formatted conversions retain their existing Rev semantics.
    void update() override
    {
        std::unique_ptr<AbstractModelObject> current(sourceObject->clone());
        current->setValueFromNode(*sourceNode);
        std::unique_ptr<RevObject> converted(current->convertTo(targetType));
        auto* value = copyValue(static_cast<const ModelObject<T>&>(*converted).getValue());
        delete this->value;
        this->value = value;
    }

protected:
    void swapParameterInternal(const RevBayesCore::DagNode*, const RevBayesCore::DagNode* replacement) override
    {
        sourceNode = replacement;
    }

private:
    static T* copyValue(const T& value)
    {
        return RevBayesCore::Cloner<T, IsDerivedFrom<T, RevBayesCore::Cloneable>::Is>::createClone(value);
    }

    const RevBayesCore::DagNode* sourceNode;
    std::shared_ptr<const AbstractModelObject> sourceObject;
    TypeSpec targetType;
};

}

#endif
