#ifndef _TIMESTAMP_TABLE_H
#define _TIMESTAMP_TABLE_H

#include <vlog/concepts.h>
#include <vlog/edbtable.h>

#include <vector>

class TimestampTable: public EDBTable {
    protected:
        const PredId_t predid;
        EDBLayer *layer;

    public:
        TimestampTable(PredId_t predid,
                EDBLayer *layer);

        uint8_t getArity() const {
            return 1;
        }

        bool areTermsEncoded() {
            return false;
        }

        void query(QSQQuery *query, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        size_t estimateCardinality(const Literal &query) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        size_t getCardinality(const Literal &query) {
            return 1;
        }

        size_t getCardinalityColumn(const Literal &query, uint8_t posColumn) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        bool isEmpty(const Literal &query, std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        EDBIterator *getSortedIterator(const Literal &query,
                const std::vector<uint8_t> &fields) {
            return getIterator(query);
        }

        bool getDictNumber(const char *text, const size_t sizeText,
                uint64_t &id) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        bool getDictText(const uint64_t id, char *text) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        bool getDictText(const uint64_t id, std::string &text) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        uint64_t getNTerms() {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        uint64_t getSize() {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        EDBIterator *getIterator(const Literal &query);

        void releaseIterator(EDBIterator *itr);

        virtual ~TimestampTable() {
        }
};

#endif
