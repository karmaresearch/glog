#include <vlog/text/stringtable.h>
#include <vlog/text/stringitr.h>

StringTable::StringTable(PredId_t predid,
        EDBLayer *layer,
        std::string fname) : predid(predid), layer(layer), fname(fname) {
    buffer1 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
    buffer2 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
}

void StringTable::query(QSQQuery *query, TupleTable *outputTable,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool StringTable::execFunction(uint64_t t1, const uint64_t t2) {
    bool found1 = layer->getDictText(t1, buffer1.get());
    bool found2 = layer->getDictText(t2, buffer2.get());
    if (!found1 || !found2) {
        return false;
    }

    if (fname == "containedIn") {
        auto resp = std::string(buffer2.get()).find(std::string(buffer1.get()));
        return resp != std::string::npos;
    } else if (fname == "equal") {
        return strcmp(buffer1.get(), buffer2.get()) == 0;
    }

    LOG(ERRORL) << "Function " << fname << " is unknown";
    throw 10;
}

size_t StringTable::getCardinality(const Literal &query) {
    auto t1 = query.getTermAtPos(0);
    auto t2 = query.getTermAtPos(1);
    if (!t1.isVariable() && !t2.isVariable()) {
        bool found = execFunction(t1.getValue(), t2.getValue());
        if (found)
            return 1;
        else
            return 0;
    }
    return 0;
}

size_t StringTable::getCardinalityColumn(const Literal &query,
        uint8_t posColumn) {
    return getCardinality(query);
}

bool StringTable::isQueryAllowed(const Literal &query) {
    auto t1 = query.getTermAtPos(0);
    auto t2 = query.getTermAtPos(1);
    return !t1.isVariable() && !t2.isVariable();
}

size_t StringTable::estimateCardinality(const Literal &query) {
    return getCardinality(query);
}

bool StringTable::isEmpty(const Literal &query,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *StringTable::getIterator(const Literal &query) {
    if (!isQueryAllowed(query)) {
        LOG(ERRORL) << "Not possible to obtain iterator for input query";
        throw 10;
    }
    auto t1 = query.getTermAtPos(0);
    auto t2 = query.getTermAtPos(1);
    bool found = execFunction(t1.getValue(), t2.getValue());
    return new StringIterator(predid, t1.getValue(), t1.getValue(), found);
}

EDBIterator *StringTable::getSortedIterator(const Literal &query,
        const std::vector<uint8_t> &fields) {
    return getIterator(query);
}

bool StringTable::getDictNumber(const char *text, const size_t sizeText,
        uint64_t &id) {
    return false;
}

bool StringTable::getDictText(const uint64_t id, char *text) {
    return false;
}

bool StringTable::getDictText(const uint64_t id, std::string &text) {
    return false;
}

uint64_t StringTable::getNTerms() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint64_t StringTable::getSize() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

void StringTable::releaseIterator(EDBIterator *itr) {
    delete itr;
}

StringTable::~StringTable() {
}
