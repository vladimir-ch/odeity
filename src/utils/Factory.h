#ifndef FACTORY_H
#define FACTORY_H

#include <map>

///  \defgroup FactoryGroup Factory

////////////////////////////////////////////////////////////////////////////////
///  \class DefaultFactoryError
///
///  \ingroup FactoryGroup
///  Manages the "Unknown Type" error in an object factory
////////////////////////////////////////////////////////////////////////////////

template <typename IdentifierType, class AbstractProduct>
struct DefaultFactoryError
{
    struct Exception : public std::exception
    {
        Exception(const IdentifierType& unknownId)
            : unknownId_(unknownId)
        { }
        virtual ~Exception() throw ()
        { }

        const char* what() const throw()
        {
            return "Unknown Type";
        }
        const IdentifierType getId()
        {
            return unknownId_;
        }
    private:
        IdentifierType unknownId_;
    };

    static AbstractProduct* OnUnknownType(const IdentifierType& id)
    {
        throw Exception(id);
    }
};

////////////////////////////////////////////////////////////////////////////////
///  \class Factory
///
///  \ingroup FactoryGroup
///  Implements a generic object factory.
///
////////////////////////////////////////////////////////////////////////////////
template
<
    class AbstractProduct,
    typename IdentifierType,
    typename ProductCreator = AbstractProduct* (*)(),
    template<typename, class> class FactoryErrorPolicy = DefaultFactoryError
>
class Factory : public FactoryErrorPolicy<IdentifierType, AbstractProduct>
{
public:
    bool registerCreator(const IdentifierType& id, ProductCreator creator)
    {
        return associations_.insert( typename AssocMap::value_type(id, creator)).second;
    }
    bool unregisterCreator(const IdentifierType& id)
    {
        return associations_.erase(id) == 1;
    }
    AbstractProduct* createObject(const IdentifierType& id)
    {
        typename AssocMap::const_iterator i = associations_.find(id);
        if (i != associations_.end())
        {
            return (i->second)();
        }
        return this->OnUnknownType(id);
    }
private:
    typedef std::map<IdentifierType, ProductCreator> AssocMap;
    AssocMap associations_;
};

#endif // FACTORY_H
