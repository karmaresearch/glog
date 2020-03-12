#ifndef _STRING_TABLE_H
#define _STRING_TABLE_H

#include <vlog/concepts.h>
#include <vlog/edbtable.h>

#include <vector>

class StringTable: public EDBTable {
    protected:
        const PredId_t predid;
        EDBLayer *layer;
        const std::string fname;

        virtual bool execFunction(const uint64_t t1, const uint64_t t2) {
            return false;
        }

        virtual bool execFunction(const uint64_t t1) {
            return false;
        }

    public:
        StringTable(PredId_t predid,
                EDBLayer *layer,
                std::string fname);

        virtual uint8_t getArity() const = 0;

        bool areTermsEncoded() {
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

        uint64_t getNTerms();

        void releaseIterator(EDBIterator *itr);

        uint64_t getSize();

        virtual ~StringTable();
};

#endif
