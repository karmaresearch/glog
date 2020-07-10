#include <vlog/clique/cliquetable.h>

CliqueTable::CliqueTable(PredId_t targetPredicate) :
    targetPredicate(targetPredicate) {
    }

void CliqueTable::query(QSQQuery *query, TupleTable *outputTable,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t CliqueTable::estimateCardinality(const Literal &query) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t CliqueTable::getCardinality(const Literal &query) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t CliqueTable::getCardinalityColumn(const Literal &query,
        uint8_t posColumn) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool CliqueTable::isEmpty(const Literal &query, std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *CliqueTable::getIterator(const Literal &query) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *CliqueTable::getSortedIterator(const Literal &query,
        const std::vector<uint8_t> &fields) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

void CliqueTable::releaseIterator(EDBIterator *itr) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool CliqueTable::getDictNumber(const char *text, const size_t sizeText,
        uint64_t &id) {
    return false;
}

bool CliqueTable::getDictText(const uint64_t id, char *text) {
    return false;
}

bool CliqueTable::getDictText(const uint64_t id, std::string &text) {
    return false;
}

uint64_t CliqueTable::getNTerms() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint64_t CliqueTable::getSize() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint8_t CliqueTable::getArity() const {
    return 2;
}
