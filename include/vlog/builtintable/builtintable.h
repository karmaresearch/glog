#ifndef _BUILTIN_TABLE_H
#define _BUILTIN_TABLE_H

#include <vlog/concepts.h>
#include <vlog/edbtable.h>

#include <vector>

class BuiltinTable: public EDBTable {
    protected:
        const PredId_t predid;
        EDBLayer *layer;
        const std::string fname;

        bool execFunction(const uint64_t t1);

        bool builtinFunction(Term_t *t, uint8_t *pos) {
            return execFunction(t[pos[0]]);
        }

    public:
        BuiltinTable(PredId_t predid,
                EDBLayer *layer,
                std::string fname);

        uint8_t getArity() const {
            return 1;
        }

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

        bool acceptQueriesWithFreeVariables() {
            return false;
        }

        BuiltinFunction getBuiltinFunction();

        virtual ~BuiltinTable();
};

#endif
