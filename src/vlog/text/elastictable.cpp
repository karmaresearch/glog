#include <vlog/text/elastictable.h>
#include <vlog/text/elasticitr.h>

#include <kognac/logs.h>

ElasticTable::ElasticTable(PredId_t predid,
        EDBLayer *layer,
        std::string baseurl,
        std::string field,
        std::string startRange) : predid(predid), layer(layer),
    baseurl(baseurl), field(field), startRange(stoi(startRange)) {
    }

void ElasticTable::query(QSQQuery *query, TupleTable *outputTable,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t ElasticTable::getCardinality(const Literal &query) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t ElasticTable::getCardinalityColumn(const Literal &query,
        uint8_t posColumn) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool ElasticTable::isQueryAllowed(const Literal &query) {
    auto t1 = query.getTermAtPos(0);
    return !t1.isVariable();
}


size_t ElasticTable::estimateCardinality(const Literal &query) {
    return getCardinality(query);
}

bool ElasticTable::isEmpty(const Literal &query,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *ElasticTable::getIterator(const Literal &query){
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *ElasticTable::getSortedIterator(const Literal &query,
        const std::vector<uint8_t> &fields) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool ElasticTable::getDictNumber(const char *text, const size_t sizeText,
        uint64_t &id) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool ElasticTable::getDictText(const uint64_t id, char *text) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool ElasticTable::getDictText(const uint64_t id, std::string &text) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint64_t ElasticTable::getNTerms() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint64_t ElasticTable::getSize() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

void ElasticTable::releaseIterator(EDBIterator *itr) {
    delete itr;
}

ElasticTable::~ElasticTable() {
}
