#ifndef _CLIQUE_TABLE_H
#define _CLIQUE_TABLE_H

#include <vlog/clique/cliqueiterator.h>

#include <vlog/edbtable.h>
#include <glog/gbgraph.h>

#include <vector>

class CliqueTable: public EDBTable {
    private:
        const PredId_t predid;
        const PredId_t targetPredicate;
        GBGraph *g;
        size_t step;

        bool recompute;
        std::map<size_t, std::vector<Term_t>> components;
        std::map<Term_t, size_t> term2component;
        size_t componentIDCounter;

        void computeConnectedComponents();

        CliqueIterator *iterator();

    public:
        CliqueTable(PredId_t predid, PredId_t targetPredicate);

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

        void setContext(GBGraph *g, size_t step);

        void clearContext();

        std::vector<std::pair<Term_t, Term_t>> checkNewIn(
                const Literal &l1,
                std::vector<uint8_t> &posInL1,
                const std::vector<std::pair<Term_t, Term_t>> &existing);
};

#endif
