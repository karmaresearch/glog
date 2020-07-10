#ifndef _CLIQUE_TABLE_H
#define _CLIQUE_TABLE_H

#include <vlog/edbtable.h>
#include <vlog/gbchase/gbgraph.h>

#include <vector>

class CliqueTable: public EDBTable {
    private:
        const PredId_t targetPredicate;
        GBGraph *g;
        size_t step;

    public:
        CliqueTable(PredId_t targetPredicate);

        void query(QSQQuery *query, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        size_t estimateCardinality(const Literal &query);

        size_t getCardinality(const Literal &query);

        size_t getCardinalityColumn(const Literal &query, uint8_t posColumn);

        bool isEmpty(const Literal &query, std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        EDBIterator *getIterator(const Literal &query);

        EDBIterator *getSortedIterator(const Literal &query,
                const std::vector<uint8_t> &fields);

        void releaseIterator(EDBIterator *itr);

        bool getDictNumber(const char *text, const size_t sizeText,
                uint64_t &id);

        bool getDictText(const uint64_t id, char *text);

        bool getDictText(const uint64_t id, std::string &text);

        uint64_t getNTerms();

        uint64_t getSize();

        uint8_t getArity() const;

        bool canChange() {
            return true;
        }

        void setContext(GBGraph *g, size_t step) {
            this->g = g;
            this->step = step;
        }
};

#endif
