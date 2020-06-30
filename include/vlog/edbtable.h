#ifndef _EDB_TABLE_H
#define _EDB_TABLE_H

#include <vlog/qsqquery.h>
#include <vlog/idxtupletable.h>
#include <vlog/builtin.h>

class Column;
class EDBIterator;
class Segment;
class TGSegment;

class EDBTable {
    public:
        virtual void join(std::vector<Term_t> &out, const Literal &l1,
                std::vector<uint8_t> &posInL1, const uint8_t joinLeftVarPos,
                const Literal &l2, const uint8_t posInL2,
                const uint8_t copyVarPosLeft);

        virtual void join(std::vector<std::pair<Term_t,Term_t>> &out,
                const Literal &l1, std::vector<uint8_t> &posInL1,
                const uint8_t joinLeftVarPos,
                const Literal &l2, const uint8_t posInL2,
                const uint8_t copyVarPosLeft1,
                const uint8_t copyVarPosLeft2);

        virtual std::vector<std::shared_ptr<Column>> checkNewIn(const Literal &l1,
                std::vector<uint8_t> &posInL1,
                const Literal &l2,
                std::vector<uint8_t> &posInL2);

        virtual std::vector<std::shared_ptr<Column>> checkNewIn(
                std::vector<std::shared_ptr<Column>> &checkValues,
                const Literal &l2,
                std::vector<uint8_t> &posInL2);

        virtual std::vector<std::pair<Term_t, Term_t>> checkNewIn(
                const Literal &l1,
                std::vector<uint8_t> &posInL1,
                const std::vector<std::pair<Term_t, Term_t>> &existing);

        virtual std::vector<Term_t> checkNewIn(
                std::shared_ptr<const TGSegment> newSeg,
                int posNew,
                const Literal &l2,
                int posInL2);

        virtual std::vector<Term_t> checkNewIn(
                const Literal &l1,
                int posInL1,
                std::shared_ptr<const TGSegment> oldSeg);

        virtual std::shared_ptr<Column> checkIn(
                const std::vector<Term_t> &values,
                const Literal &l2,
                uint8_t posInL2,
                size_t &sizeOutput);

        //execute the query on the knowledge base
        virtual void query(QSQQuery *query, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter) = 0;

        virtual size_t estimateCardinality(const Literal &query) = 0;

        virtual size_t getCardinality(const Literal &query) = 0;

        virtual size_t getCardinalityColumn(const Literal &query, uint8_t posColumn) = 0;

        virtual bool isEmpty(const Literal &query, std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter) = 0;

        virtual EDBIterator *getIterator(const Literal &query) = 0;

        virtual EDBIterator *getSortedIterator(const Literal &query,
                const std::vector<uint8_t> &fields) = 0;

        virtual void releaseIterator(EDBIterator *itr) = 0;

        virtual bool getDictNumber(const char *text, const size_t sizeText,
                uint64_t &id) = 0;

        virtual bool getDictText(const uint64_t id, char *text) = 0;

        virtual bool getDictText(const uint64_t id, std::string &text) = 0;

        virtual uint64_t getNTerms() = 0;

        virtual uint64_t getSize() = 0;

        virtual bool areTermsEncoded() {
            return true;
        }

        virtual uint8_t getArity() const = 0;

        virtual bool useSegments() {
            return false;
        }

        virtual bool expensiveLayer() {
            return false;
        }

        virtual std::shared_ptr<const Segment> getSegment() {
            return std::shared_ptr<const Segment>();
        }

        virtual bool isQueryAllowed(const Literal &query) {
            return true;
        }

        virtual bool acceptQueriesWithFreeVariables() {
            return true;
        }

        virtual BuiltinFunction getBuiltinFunction() {
            LOG(ERRORL) << "This function should have been implemented"
                " in a subclass";
            throw 10;
        }
};


#endif
