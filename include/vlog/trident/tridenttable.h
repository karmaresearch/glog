#ifndef _TRIDENT_TABLE_H
#define _TRIDENT_TABLE_H

#include <vlog/trident/tridentiterator.h>
#include <vlog/column.h>
#include <vlog/edbtable.h>


#include <trident/kb/kb.h>
#include <trident/kb/querier.h>
#include <kognac/factory.h>

class SeqColumnWriter : public SequenceWriter {
    private:
        ColumnWriter *writer;
    public:
        SeqColumnWriter(ColumnWriter *w) : writer(w) {
        }

        void add(const uint64_t v) {
            writer->add(v);
        }
};

class TridentTable: public EDBTable {

    private:

        KB *kb;
        Querier *q;
        DictMgmt *dict;
        Factory<TridentIterator> kbItrFactory;
        std::mutex mutex;
        bool multithreaded;
        EDBLayer *layer;

        TridentIterator *getTridentIter();

        std::vector<std::shared_ptr<Column>> performAntiJoin(const Literal &l1,
                std::vector<uint8_t> &pos1, const Literal &l2,
                std::vector<uint8_t> &pos2,
                bool stopAfterFirst = false);

        std::vector<std::pair<Term_t, Term_t>> performAntiJoin(const Literal &l1,
                std::vector<uint8_t> &pos1,
                const std::vector<std::pair<Term_t, Term_t>> &existing);

        std::vector<std::shared_ptr<Column>> performAntiJoin(
                std::vector<std::shared_ptr<Column>>
                &valuesToCheck,
                const Literal &l,
                std::vector<uint8_t> &pos);

        void getQueryFromEDBRelation0(QSQQuery *query, TupleTable *outputTable);

        void getQueryFromEDBRelation12(QSQQuery *query, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        void getQueryFromEDBRelation12(VTerm s, VTerm p, VTerm o, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<std::pair<uint64_t, uint64_t>> *pairs,
                std::vector<int> &posVarsToReturn,
                std::vector<std::pair<int, int>> &joins,
                std::vector<std::vector<int>> &posToCopy);

        void getQueryFromEDBRelation3(QSQQuery *query, TupleTable *outputTable,
                std::vector<Term_t> *valuesToFilter);

        void join(std::vector<Term_t> &out1,
                std::vector<std::pair<Term_t,Term_t>> &out2,
                const Literal &l1,
                std::vector<uint8_t> &posInL1, const uint8_t joinLeftVarPos,
                const Literal &l2, const uint8_t posInL2,
                const std::vector<uint8_t> copyVarPosLeft);

    public:
        TridentTable(std::string kbDir, bool multithreaded, EDBLayer *layer) : layer(layer) {
            KBConfig config;
            kb = new KB(kbDir.c_str(), true, false, true, config);
            q = kb->query();
            dict = kb->getDictMgmt();
            this->multithreaded = multithreaded;
        }

        void join(std::vector<Term_t> &out, const Literal &l1,
                std::vector<uint8_t> &posInL1, const uint8_t joinLeftVarPos,
                const Literal &l2, const uint8_t posInL2,
                const uint8_t copyVarPosLeft) {
            std::vector<std::pair<Term_t, Term_t>> unused;
            std::vector<uint8_t> pos;
            pos.push_back(copyVarPosLeft);
            join(out, unused, l1, posInL1, joinLeftVarPos,
                    l2, posInL2, pos);
        }

        void join(std::vector<std::pair<Term_t,Term_t>> &out,
                const Literal &l1, std::vector<uint8_t> &posInL1,
                const uint8_t joinLeftVarPos,
                const Literal &l2, const uint8_t posInL2,
                const uint8_t copyVarPosLeft1,
                const uint8_t copyVarPosLeft2) {
            std::vector<Term_t> unused;
            std::vector<uint8_t> pos;
            pos.push_back(copyVarPosLeft1);
            pos.push_back(copyVarPosLeft2);
            join(unused, out, l1, posInL1, joinLeftVarPos,
                    l2, posInL2, pos);
        }

        std::vector<std::shared_ptr<Column>> checkNewIn(
                const Literal &l1,
                std::vector<uint8_t> &posInL1,
                const Literal &l2,
                std::vector<uint8_t> &posInL2,
                bool stopAfterFirst = false);

        std::vector<std::shared_ptr<Column>> checkNewIn(
                std::vector <
                std::shared_ptr<Column >> &checkValues,
                const Literal &l2,
                std::vector<uint8_t> &posInL2);

        std::vector<std::pair<Term_t, Term_t>> checkNewIn(
                const Literal &l1,
                std::vector<uint8_t> &posInL1,
                const std::vector<std::pair<Term_t, Term_t>> &existing);

        std::vector<Term_t> checkNewIn(
                std::shared_ptr<const TGSegment> newSeg,
                int posNew,
                const Literal &l2,
                int posInL2);

        std::vector<Term_t> checkNewIn(
                const Literal &l1,
                int posInL1,
                std::shared_ptr<const TGSegment> oldSeg);

        std::vector<std::pair<Term_t,Term_t>> checkNewIn(
                std::shared_ptr<const TGSegment> newSeg,
                int posNew1,
                int posNew2,
                const Literal &l2,
                int posInL2_1,
                int posInL2_2);

        std::vector<std::pair<Term_t, Term_t>> checkNewIn(
                const std::vector<std::pair<Term_t, Term_t>> &terms,
                const Literal &l2,
                int posInL2_1,
                int posInL2_2);

        std::shared_ptr<Column> checkIn(
                const std::vector<Term_t> &values,
                const Literal &l2,
                uint8_t posInL2,
                size_t &sizeOutput);

        //execute the query on the knowledge base
        void query(QSQQuery *query, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        Querier *getQuerier() {
            return q;
        }

        KB *getKB() {
            return kb;
        }

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

        uint8_t getArity() const {
            // Moved this manifest constant from edb importer to here       RFHH
            return 3;
        }

        uint64_t getNTerms();

        uint64_t getSize();

        void releaseIterator(EDBIterator *itr);

        ~TridentTable() {
            delete q;
            delete kb;
        }

};

#endif
