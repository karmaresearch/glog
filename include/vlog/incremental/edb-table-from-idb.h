#ifndef INCREMENTAL__EDB_TABLE_FROM_IDB_H__
#define INCREMENTAL__EDB_TABLE_FROM_IDB_H__

/**
 * In the overdelete and rederive steps DRed algorithm, we require the values
 * of the previous materialization stages, but they are immutable. So we treat
 * them as if they are EDB tables.
 */

#include <vlog/fcinttable.h>
#include <vlog/seminaiver.h>
#include <vlog/inmemory/inmemorytable.h>

class EDBonIDBIterator : public EDBIterator {
protected:
    const Literal &query;
    const std::shared_ptr<SemiNaiver> prevSemiNaiver;
    const PredId_t predid;
    FCIterator idbItr;
    std::shared_ptr<const FCInternalTable> idbInternalTable;
    FCInternalTableItr *idbInternalItr;

public:
    EDBonIDBIterator(const Literal &query, const std::shared_ptr<SemiNaiver> SN) :
            query(query), prevSemiNaiver(SN),
            predid(query.getPredicate().getId()), idbInternalItr(NULL) {
        init();
    }

    virtual void init() {
        idbItr = prevSemiNaiver->getTable(query, 0,
                                          prevSemiNaiver->getCurrentIteration());
        if (! idbItr.isEmpty()) {
            idbInternalTable = idbItr.getCurrentTable();
            idbInternalItr = idbInternalTable->getIterator();
        }
    }

    virtual bool hasNext() {
        LOG(DEBUGL) << "         hasNext() on " << query.tostring();
        if (idbInternalItr == NULL) {
            return false;
        }
        if (idbInternalItr->hasNext()) {
            return true;
        }

        idbInternalTable->releaseIterator(idbInternalItr);
        idbItr.moveNextCount();
        if (idbItr.isEmpty()) {
            idbInternalItr = NULL;
            return false;
        }

        idbInternalTable = idbItr.getCurrentTable();
        idbInternalItr = idbInternalTable->getIterator();

        return idbInternalItr->hasNext();
    }

    virtual void next() {
        LOG(DEBUGL) << "         next() on " << query.tostring();
        if (! idbInternalItr->hasNext()) {
            LOG(ERRORL) << "next() on a table that hasNext() = false";
            throw 10;               // follow convention
        }
        return idbInternalItr->next();
    }

    virtual Term_t getElementAt(const uint8_t p) {
        LOG(DEBUGL) << "get element[" << (int)p << "] of " << query.tostring() << " = " << idbInternalItr->getCurrentValue(p);
        return idbInternalItr->getCurrentValue(p);
    }

    virtual PredId_t getPredicateID() {
        return predid;
    }

    virtual void moveTo(const uint8_t field, const Term_t t) {
        throw 10;
    }

    virtual void skipDuplicatedFirstColumn() {
        LOG(ERRORL) << "FIXME: implement " << typeid(*this).name() << "::" <<
            __func__;
    }

    virtual void clear() {
        LOG(ERRORL) << "FIXME: implement " << typeid(*this).name() << "::" <<
            __func__;
    }

    virtual const char *getUnderlyingArray(uint8_t column) {
        return NULL;
    }

    /*
    No need to override
    virtual std::pair<uint8_t, std::pair<uint8_t, uint8_t>> getSizeElemUnderlyingArray(uint8_t column) {
        return std::make_pair(0, make_pair(0, 0));
    }
    */

    virtual ~EDBonIDBIterator() {}
};

class EDBonIDBSortedIterator : public EDBIterator {
protected:
    const Literal &query;
    const std::vector<uint8_t> &fields;
    const std::shared_ptr<SemiNaiver> prevSemiNaiver;
    EDBLayer *layer;
    const PredId_t predid;
    InmemoryTable *inmemoryTable;
    EDBIterator *itr;

public:
    EDBonIDBSortedIterator(const Literal &query,
                           const std::vector<uint8_t> &fields,
                           const std::shared_ptr<SemiNaiver> SN,
                           EDBLayer *layer) :
            query(query), fields(fields), prevSemiNaiver(SN), layer(layer),
            predid(query.getPredicate().getId()) {
        LOG(DEBUGL) << "EDBonIDBSortedIterator, query " << query.tostring() <<
            ", fields.size() " << fields.size();
        LOG(INFOL) << "FIXME: if fields implies a no-op, dispatch to non-sorted Iterator";
        EDBonIDBIterator between(query, prevSemiNaiver);
        // Collect matching data. Will be stored in an InmemoryTable.
        std::vector<std::vector<Term_t>> rows;
        std::vector<Term_t> term(query.getTupleSize());
        // need to store all variables, then afterwards sort by fields
        std::vector<bool> isVariable(query.getTupleSize());
        VTuple tuple = query.getTuple();
        for (size_t i = 0; i < query.getTupleSize(); ++i) {
            VTerm vt = tuple.get(i);
            isVariable[i] = vt.isVariable();
            if (! isVariable[i]) {
                // fill out the constants
                term[i] = vt.getValue();
            }
        }
        while (between.hasNext()) {
            between.next();
            size_t nVars = 0;
            for (size_t i = 0; i < query.getTupleSize(); ++i) {
                if (isVariable[i]) {
                    term[i] = between.getElementAt(nVars);
                    ++nVars;
                }
            }
            rows.push_back(term);
        }

        if (rows.size() == 0) {
            inmemoryTable = NULL;
        } else {
            inmemoryTable = new InmemoryTable(predid, rows, layer);
            itr = inmemoryTable->getSortedIterator(query, fields);
        }
    }

    virtual bool hasNext() {
        LOG(DEBUGL) << "         hasNext() on " << query.tostring();
        if (inmemoryTable == NULL) {
            return false;
        }

        return itr->hasNext();
    }

    virtual void next() {
        LOG(DEBUGL) << "         next() on " << query.tostring();
        return itr->next();
    }

    virtual Term_t getElementAt(const uint8_t p) {
        LOG(DEBUGL) << "get element[" << (int)p << "] of " << query.tostring() << " = " << itr->getElementAt(p);
        return itr->getElementAt(p);
    }

    virtual PredId_t getPredicateID() {
        return predid;
    }

    virtual void moveTo(const uint8_t field, const Term_t t) {
        throw 10;
    }

    virtual void skipDuplicatedFirstColumn() {
        LOG(ERRORL) << "FIXME: implement " << typeid(*this).name() << "::" <<
            __func__;
    }

    virtual void clear() {
        LOG(ERRORL) << "FIXME: implement " << typeid(*this).name() << "::" <<
            __func__;
    }

    virtual const char *getUnderlyingArray(uint8_t column) {
        return NULL;
    }

    /*
    No need to override
    virtual std::pair<uint8_t, std::pair<uint8_t, uint8_t>> getSizeElemUnderlyingArray(uint8_t column) {
        return std::make_pair(0, make_pair(0, 0));
    }
    */

    virtual ~EDBonIDBSortedIterator() {
        LOG(ERRORL) << "Need to cache inmemoryTable";
        if (inmemoryTable != NULL) {
            inmemoryTable->releaseIterator(itr);
            delete inmemoryTable;
        }
    }
};

/**
 * Transfer an IDB predicate of the previous SemiNaiver to be
 * a new 'EDB' predicate.
 */
class EDBonIDBTable : public EDBTable {
private:
    PredId_t    predid;
    const std::shared_ptr<SemiNaiver> prevSemiNaiver;
    const FCTable *idbTable;
    size_t cardinality;

public:
    EDBLayer *layer;
    EDBonIDBTable(PredId_t predid, EDBLayer *layer,
                  const std::shared_ptr<SemiNaiver> prevSN) :
            predid(predid), layer(layer), prevSemiNaiver(prevSN),
            idbTable(prevSN->getTable(predid, prevSN->getCurrentIteration())) {
        cardinality = prevSN->getProgram()->getPredicateCard(predid);
        LOG(DEBUGL) << "prevSN iteration " << prevSN->getCurrentIteration() <<
            " card. " << cardinality << " size " << idbTable->getNAllRows();
    }

    //execute the query on the knowledge base
    virtual void query(QSQQuery *query, TupleTable *outputTable,
                       std::vector<uint8_t> *posToFilter,
                       std::vector<Term_t> *valuesToFilter) {
        LOG(ERRORL) << "FIXME: implement " << __func__;
    }

    virtual size_t estimateCardinality(const Literal &query) {
        return getCardinality(query);
    }

    virtual size_t getCardinality(const Literal &query) {
        return cardinality;
    }

    virtual size_t getCardinalityColumn(const Literal &query, uint8_t posColumn) {
        LOG(ERRORL) << "FIXME: implement " << __func__;
    }

    virtual bool isEmpty(const Literal &query, std::vector<uint8_t> *posToFilter,
                         std::vector<Term_t> *valuesToFilter) {
       if (posToFilter == NULL) {
           return getCardinality(query) == 0;
       }

       LOG(ERRORL) << "Not implemented:" << __func__;
       throw 10;
    }

    virtual EDBIterator *getIterator(const Literal &q) {
        LOG(DEBUGL) << "Get iterator for query " << q.tostring(NULL, layer);
        return new EDBonIDBIterator(q, prevSemiNaiver);
    }

    virtual EDBIterator *getSortedIterator(const Literal &query,
                                           const std::vector<uint8_t> &fields) {
        LOG(INFOL) << "Get SortedIterator for query " <<
            query.tostring(NULL, layer) << " fields " << fields.size();
        LOG(INFOL) << "FIXME: implement cache for sorted query";
        return new EDBonIDBSortedIterator(query, fields, prevSemiNaiver, layer);
    }

    virtual void releaseIterator(EDBIterator *itr) {
        delete itr;
    }

    virtual bool getDictNumber(const char *text, const size_t sizeText,
                               uint64_t &id) {
        return false;
    }

    virtual bool getDictText(const uint64_t id, char *text) {
        return false;
    }

    virtual bool getDictText(const uint64_t id, std::string &text) {
        return false;
    }

    virtual uint64_t getNTerms() {
        return 0;
    }

    virtual uint8_t getArity() const {
        return idbTable->getSizeRow();
    }

    virtual uint64_t getSize() {
        return idbTable->getNAllRows();
    }

    virtual PredId_t getPredicateID() const {
        return predid;
    }

    std::ostream &dump(std::ostream &os) {
        std::string name = layer->getPredName(predid);
        os << (void *)this << " Table " << name << std::endl;
        const uint8_t sizeRow = getArity();
        Predicate pred(predid, 0, EDB, sizeRow);
        VTuple t = VTuple(sizeRow);
        for (uint8_t i = 0; i < t.getSize(); ++i) {
            t.set(VTerm(i + 1, 0), i);
        }
        Literal lit(pred, t);
        os << "Literal@" << int(sizeRow) << " " << lit.tostring() << std::endl;
        EDBIterator *itr = getIterator(lit);
        while (itr->hasNext()) {
            itr->next();
            for (uint8_t m = 0; m < sizeRow; ++m) {
                os << "\t";
                uint64_t v = itr->getElementAt(m);
                std::string buffer;
                if (getDictText(v, buffer)) {
                    os << buffer;
                } else {
                    os << to_string(v);
                }
            }
            os << std::endl;
        }

        releaseIterator(itr);

        return os;
    }
};

/*
std::ostream &operator<<(ostream &os, const EDBonIDBTable &table) {
    return os << table.dump();
}
*/

#endif  // INCREMENTAL__EDB_TABLE_FROM_IDB_H__
