#include <vlog/builtintable/builtintable.h>
#include <vlog/chasemgmt.h>
#include <vlog/text/stringitr.h>

BuiltinTable::BuiltinTable(PredId_t predid,
        EDBLayer *layer,
        std::string fname) : predid(predid), layer(layer), fname(fname) {
}

void BuiltinTable::query(QSQQuery *query, TupleTable *outputTable,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t BuiltinTable::getCardinality(const Literal &query) {
    auto t1 = query.getTermAtPos(0);
    if (!t1.isVariable()) {
        bool found = execFunction(t1.getValue());
        if (found)
            return 1;
        else
            return 0;
    }
    return 0;
}

size_t BuiltinTable::getCardinalityColumn(const Literal &query,
        uint8_t posColumn) {
    return getCardinality(query);
}

bool BuiltinTable::isQueryAllowed(const Literal &query) {
    for(int i = 0; i < getArity(); ++i) {
        auto t = query.getTermAtPos(i);
        if (t.isVariable())
            return false;
    }
    return true;
}

size_t BuiltinTable::estimateCardinality(const Literal &query) {
    return getCardinality(query);
}

bool BuiltinTable::isEmpty(const Literal &query,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *BuiltinTable::getIterator(const Literal &query) {
    if (!isQueryAllowed(query)) {
        LOG(ERRORL) << "Not possible to obtain iterator for input query";
        throw 10;
    }
    auto t1 = query.getTermAtPos(0);
    bool found = execFunction(t1.getValue());
    return new StringIterator(predid, t1.getValue(), found);
}

EDBIterator *BuiltinTable::getSortedIterator(const Literal &query,
        const std::vector<uint8_t> &fields) {
    return getIterator(query);
}

bool BuiltinTable::getDictNumber(const char *text, const size_t sizeText,
        uint64_t &id) {
    return false;
}

bool BuiltinTable::getDictText(const uint64_t id, char *text) {
    return false;
}

bool BuiltinTable::getDictText(const uint64_t id, std::string &text) {
    return false;
}

uint64_t BuiltinTable::getNTerms() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint64_t BuiltinTable::getSize() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

void BuiltinTable::releaseIterator(EDBIterator *itr) {
    delete itr;
}

bool BuiltinTable::execFunction(const uint64_t t1)
{
    if (this->fname=="IS_NULL")
    {
        return IS_NULLVALUE(t1);
    }
    return false;
}

BuiltinFunction BuiltinTable::getBuiltinFunction() {
    BuiltinFunction fn;
    fn.fn = std::bind(&BuiltinTable::builtinFunction,
            this,
            std::placeholders::_1, std::placeholders::_2);
    return fn;
}

BuiltinTable::~BuiltinTable() {
}
