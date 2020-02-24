#ifndef _ELASTIC_TABLE_H
#define _ELASTIC_TABLE_H

#include <vlog/concepts.h>
#include <vlog/edbtable.h>

#include <vector>

class ElasticTable: public EDBTable {
    private:
        const PredId_t predid;
        EDBLayer *layer;
        const std::string baserel;
        const std::string basehost;
        const uint32_t baseport;
        const std::string basepath;
        const int64_t startRange;
        size_t nterms;
        std::shared_ptr<EDBTable> dictTable;

    public:
        ElasticTable(PredId_t predid,
                EDBLayer *layer,
                std::string baserel,
                std::string basehost,
                std::string baseport,
                std::string basepath,
                std::string startRange);

        uint8_t getArity() const {
            return 2;
        }

        bool areTermsEncoded() {
            return true;
        }

        bool expensiveLayer() {
            return true;
        }

        void query(QSQQuery *query, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        bool isQueryAllowed(const Literal &query);

        size_t estimateCardinality(const Literal &query);

        size_t getCardinality(const Literal &query);

        size_t getCardinalityColumn(const Literal &query, uint8_t posColumn);

        bool isEmpty(const Literal &query, std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        EDBIterator *getIterator(const Literal &query);

        EDBIterator *getSortedIterator(const Literal &query,
                const std::vector<uint8_t> &fields);

        bool getDictNumber(const char *text, const size_t sizeText,
                uint64_t &id);

        bool getDictText(const uint64_t id, char *text);

        bool getDictText(const uint64_t id, std::string &text);

        uint64_t getStartOffset() {
            return startRange;
        }

        uint64_t getNTerms();

        void releaseIterator(EDBIterator *itr);

        uint64_t getSize();

        ~ElasticTable();
};

#endif
