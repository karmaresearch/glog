#ifndef EDB_LAYER_H
#define EDB_LAYER_H

#include <vlog/consts.h>
#include <vlog/concepts.h>
#include <vlog/qsqquery.h>
#include <vlog/support.h>
#include <vlog/idxtupletable.h>

#include <vlog/edbtable.h>
#include <vlog/edbiterator.h>
#include <vlog/edbconf.h>

#include <vlog/incremental/removal.h>

#include <kognac/factory.h>

#include <vector>
#include <map>

//Datatype is set in the most significant three bits
#define IS_NUMBER(x) ((x) >> 61)
#define IS_UINT(x) ((x >> 61) == 1)
#define IS_FLOAT32(x) ((x >> 61) == 2)
#define GET_UINT(x) (x & 0x1FFFFFFFFFFFFFFFul)
#define GET_FLOAT32(x) (*(float*)&x)
#define FLOAT32_MASK(x) ( *((uint32_t*)&x) | 0x4000000000000000ul)
#define UINT_MASK(x) ( *((uint32_t*)&x) | 0x2000000000000000ul)

class Column;
class SemiNaiver;       // Why cannot I break the software hierarchy? RFHH
class EDBFCInternalTable;
class TGSegment;
class GBGraph;

using RemoveLiteralOf = std::unordered_map<PredId_t, const EDBRemoveLiterals *>;

using NamedSemiNaiver = std::unordered_map<std::string,
      std::shared_ptr<SemiNaiver>>;

class EDBMemIterator : public EDBIterator {
    private:
        uint8_t nfields = 0;
        bool isFirst = false, hasFirst = false;
        bool equalFields = false, isNextCheck = false, isNext = false;
        bool ignoreSecondColumn = false;
        bool isIgnoreAllowed = true;
        PredId_t predid;

        std::vector<Term_t>::iterator oneColumn;
        std::vector<Term_t>::iterator endOneColumn;

        std::vector<std::pair<Term_t, Term_t>>::iterator pointerEqualFieldsNext;
        std::vector<std::pair<Term_t, Term_t>>::iterator twoColumns;
        std::vector<std::pair<Term_t, Term_t>>::iterator endTwoColumns;

    public:
        EDBMemIterator() {}

        void init1(PredId_t id, std::vector<Term_t>*, const bool c1,
                const Term_t vc1);

        void init2(PredId_t id, const bool defaultSorting,
                std::vector<std::pair<Term_t, Term_t>>*, const bool c1,
                const Term_t vc1, const bool c2, const Term_t vc2,
                const bool equalFields);

        VLIBEXP void skipDuplicatedFirstColumn();

        VLIBEXP bool hasNext();

        VLIBEXP void next();

        PredId_t getPredicateID() {
            return predid;
        }

        void reset();

        void mark();

        void moveTo(const uint8_t fieldId, const Term_t t) {}

        VLIBEXP Term_t getElementAt(const uint8_t p);

        void clear() {}

        ~EDBMemIterator() {}
};

class EmptyEDBIterator final : public EDBIterator {
    PredId_t predid;

    public:
    EmptyEDBIterator(PredId_t id) {
        predid = id;
    }

    VLIBEXP void init1(PredId_t id, std::vector<Term_t>*, const bool c1,
            const Term_t vc1) {
        predid = id;
    }

    VLIBEXP void init2(PredId_t id, const bool defaultSorting,
            std::vector<std::pair<Term_t, Term_t>>*, const bool c1,
            const Term_t vc1, const bool c2, const Term_t vc2,
            const bool equalFields) {
        predid = id;
    }

    VLIBEXP void skipDuplicatedFirstColumn() { }

    VLIBEXP bool hasNext() {
        return false;
    }

    VLIBEXP void next() {
        throw 10;
    }

    PredId_t getPredicateID() {
        return predid;
    }

    void moveTo(const uint8_t fieldId, const Term_t t) {}

    VLIBEXP Term_t getElementAt(const uint8_t p) {
        throw 10;
    }

    void clear() {}

    ~EmptyEDBIterator() {}
};

class EDBLayer {
    private:

        struct EDBInfoTable {
            PredId_t id;
            uint8_t arity;
            std::string type;
            std::shared_ptr<EDBTable> manager;
        };

        const EDBConf &conf;
        const bool multithreaded;
        const std::string edbconfpath;
        bool loadAllData;

        GBGraph *context_gbGraph;
        size_t context_step;

        std::shared_ptr<Dictionary> predDictionary; //std::string, Term_t
        std::map<PredId_t, EDBInfoTable> dbPredicates;

        Factory<EDBMemIterator> memItrFactory;
        std::vector<IndexedTupleTable *>tmpRelations;

        std::vector<std::shared_ptr<EDBTable>> edbTablesWithDict;
        std::shared_ptr<Dictionary> termsDictionary;//std::string, Term_t
        std::string rootPath;

        VLIBEXP void addTridentTable(const EDBConf::Table &tableConf,
                bool multithreaded);

        void addCliqueTable(const EDBConf::Table &tableConf,
                PredId_t predId = 0,
                bool usePredId = false);

        VLIBEXP void addTopKTable(const EDBConf::Table &tableConf);

        VLIBEXP void addEmbTable(const EDBConf::Table &tableConf);

        void addElasticTable(const EDBConf::Table &tableConf);

        void addStringTable(bool isUnary, const EDBConf::Table &tableConf);

        void addBuiltinTable(const EDBConf::Table &tableConf);

        void addTimestampTable(const EDBConf::Table &tableConf);

        void addInfRoundTable(const EDBConf::Table &tableConf);

#ifdef MYSQL
        VLIBEXP void addMySQLTable(const EDBConf::Table &tableConf);
#endif
#ifdef ODBC
        VLIBEXP void addODBCTable(const EDBConf::Table &tableConf);
#endif
#ifdef MAPI
        VLIBEXP void addMAPITable(const EDBConf::Table &tableConf);
#endif
#ifdef MDLITE
        VLIBEXP void addMDLiteTable(const EDBConf::Table &tableConf);
#endif
        VLIBEXP void addInmemoryTable(const EDBConf::Table &tableConf,
                std::string edbconfpath);
        VLIBEXP void addSparqlTable(const EDBConf::Table &tableConf);

        VLIBEXP void addEDBonIDBTable(const EDBConf::Table &tableConf);
        VLIBEXP void addEDBimporter(const EDBConf::Table &tableConf);

        // literals to be removed during iteration
        RemoveLiteralOf removals;

        // For incremental reasoning: EDBonIDB must know which SemiNaiver(s) to
        // query
        NamedSemiNaiver prevSemiNaiver;

        // need to import the mapping predid -> Predicate from prevSemiNaiver
        VLIBEXP void handlePrevSemiNaiver();

        std::string name;

        void addTable(const EDBConf::Table &table, bool multithreaded,
                std::string edbconfpath, PredId_t predId = 0,
                bool usePredId = false);

    public:
        EDBLayer(EDBLayer &db, bool copyTables = false);

        EDBLayer(const EDBConf &conf, bool multithreaded,
                const NamedSemiNaiver &prevSemiNaiver,
                bool loadAllData = true) :
            conf(conf), prevSemiNaiver(prevSemiNaiver),
            loadAllData(loadAllData), multithreaded(multithreaded),
            edbconfpath(conf.getConfigFilePath()), context_gbGraph(NULL),
            context_step(0)
    {

        const std::vector<EDBConf::Table> tables = conf.getTables();
        rootPath = conf.getRootPath();
        std::string edbconfpath = conf.getConfigFilePath();

        predDictionary = std::shared_ptr<Dictionary>(new Dictionary());

        if (prevSemiNaiver.size() != 0) {
            handlePrevSemiNaiver();
        }

        for (const auto &table : tables) {
            addTable(table, multithreaded, edbconfpath);
        }
    }

        EDBLayer(const EDBConf &conf, bool multithreaded,
                bool loadAllData = true) :
            EDBLayer(conf, multithreaded, NamedSemiNaiver(), loadAllData)
    {
    }

        std::vector<PredId_t> getAllEDBPredicates();

        std::vector<PredId_t> getAllPredicateIDs() const;

        std::string getPredType(PredId_t id) const;

        VLIBEXP uint64_t getPredSize(PredId_t id) const;

        std::string getPredName(PredId_t id) const;

        uint8_t getPredArity(PredId_t id) const;

        PredId_t getPredID(const std::string &name) const;

        void setPredArity(PredId_t id, uint8_t arity);

        void addTmpRelation(Predicate &pred, IndexedTupleTable *table);

        void addEDBPredicate(std::string name,
                std::string type,
                std::vector<std::string> args,
                PredId_t id);

        PredId_t addEDBPredicate(std::string predName);

        bool isTmpRelationEmpty(Predicate &pred) {
            if (pred.getId() >= tmpRelations.size()) {
                return false;
            }
            return tmpRelations[pred.getId()] == NULL ||
                tmpRelations[pred.getId()]->getNTuples() == 0;
        }

        const Dictionary &getPredDictionary() {
            assert(predDictionary.get() != NULL);
            return *(predDictionary.get());
        }

        std::unordered_map<PredId_t, uint8_t> getPredicateCardUnorderedMap () {
            std::unordered_map< PredId_t, uint8_t> ret;
            for (auto item : dbPredicates)
                ret.insert(std::make_pair(item.first, item.second.arity));
            return ret;
        }

        VLIBEXP bool doesPredExists(PredId_t id) const;

        PredId_t getFirstEDBPredicate() {
            if (!dbPredicates.empty()) {
                auto p = dbPredicates.begin();
                return p->first;
            } else {
                LOG(ERRORL) << "There is no EDB Predicate!";
                throw 10;
            }
        }

        bool checkValueInTmpRelation(const uint8_t relId,
                const uint8_t posInRelation,
                const Term_t value) const;

        size_t getSizeTmpRelation(Predicate &pred) {
            if (pred.getId() >= tmpRelations.size()) {
                return 0;
            }
            return tmpRelations[pred.getId()]->getNTuples();
        }

        bool supportsCheckIn(const Literal &l);

        void join(std::vector<Term_t> &out, const Literal &l1,
                std::vector<uint8_t> &posInL1, const uint8_t joinLeftVarPos,
                const Literal &l2, const uint8_t posInL2,
                const uint8_t copyVarPosLeft);

        void join(std::vector<std::pair<Term_t,Term_t>> &out,
                const Literal &l1, std::vector<uint8_t> &posInL1,
                const uint8_t joinLeftVarPos,
                const Literal &l2, const uint8_t posInL2,
                const uint8_t copyVarPosLeft1,
                const uint8_t copyVarPosLeft2);

        std::vector<std::shared_ptr<Column>> checkNewIn(const Literal &l1,
                std::vector<uint8_t> &posInL1,
                const Literal &l2,
                std::vector<uint8_t> &posInL2,
                bool stopAfterFirst = false);

        // Note: posInL2 contains positions in the literal.
        std::vector<std::shared_ptr<Column>> checkNewIn(
                std::vector <
                std::shared_ptr<Column >> &checkValues,
                const Literal &l2,
                std::vector<uint8_t> &posInL2);

        std::vector<Term_t> checkNewIn(
                std::shared_ptr<const TGSegment> newSeg,
                int posNew,
                const Literal &l2,
                int posInL2);

        std::vector<Term_t> checkNewIn(
                const Literal &l1,
                int posInL1,
                std::shared_ptr<const TGSegment> oldSeg);

        std::vector<std::pair<Term_t, Term_t>> checkNewIn(const Literal &l1,
                std::vector<uint8_t> &posInL1,
                const std::vector<std::pair<Term_t, Term_t>> &existing);

        std::vector<std::pair<Term_t,Term_t>> checkNewIn(
                std::shared_ptr<const TGSegment> newSeg,
                int posNew1,
                int posNew2,
                const Literal &l2,
                int posInL2_1,
                int posInL2_2);

        std::vector<std::pair<Term_t,Term_t>> checkNewIn(
                const std::vector<std::pair<Term_t, Term_t>> &terms,
                const Literal &l2,
                int posInL2_1,
                int posInL2_2);

        std::shared_ptr<Column> checkIn(
                const std::vector<Term_t> &values,
                const Literal &l2,
                uint8_t posInL2,
                size_t &sizeOutput);

        VLIBEXP void query(QSQQuery *query, TupleTable *outputTable,
                std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        EDBIterator *getIterator(const Literal &query);

        // Note: fields only counts variables in the query.
        EDBIterator *getSortedIterator(const Literal &query,
                const std::vector<uint8_t> &fields);

        // Note: posToFilter contains positions in the literal.
        bool isEmpty(const Literal &query, std::vector<uint8_t> *posToFilter,
                std::vector<Term_t> *valuesToFilter);

        size_t estimateCardinality(const Literal &query);

        size_t getCardinality(const Literal &query);

        // Note: posColumn contains a position in the literal.
        size_t getCardinalityColumn(const Literal &query,
                uint8_t posColumn);

        VLIBEXP bool getDictNumber(const char *text,
                const size_t sizeText, uint64_t &id) const;

        VLIBEXP bool getOrAddDictNumber(const char *text,
                const size_t sizeText, uint64_t &id,
                bool detectDatatype = false);

        VLIBEXP bool getDictText(const uint64_t id, char *text) const;

        VLIBEXP bool getDictText(const uint64_t id, std::string &out) const;

        VLIBEXP std::string getDictText(const uint64_t id) const;

        Predicate getDBPredicate(int idx) const;

        std::shared_ptr<EDBTable> getEDBTable(PredId_t id) const {
            if (dbPredicates.count(id)) {
                return dbPredicates.find(id)->second.manager;
            } else {
                return std::shared_ptr<EDBTable>();
            }
        }

        std::string getTypeEDBPredicate(PredId_t id) const {
            if (dbPredicates.count(id)) {
                return dbPredicates.find(id)->second.type;
            } else {
                return "";
            }
        }

        VLIBEXP uint64_t getNTerms() const;

        VLIBEXP uint64_t getNPredicates() const;

        bool expensiveEDBPredicate(PredId_t id) {
            if (dbPredicates.count(id)) {
                return dbPredicates.find(id)->second.manager->expensiveLayer();
            }
            return false;
        }

        void releaseIterator(EDBIterator *itr);

        VLIBEXP void addRemoveLiterals(const RemoveLiteralOf &rm) {
            removals.insert(rm.begin(), rm.end());
        }

        VLIBEXP bool hasRemoveLiterals(PredId_t pred) const {
            if (removals.empty()) {
                return false;
            }
            if (removals.find(pred) == removals.end()) {
                return false;
            }
            if (removals.at(pred)->size() == 0) {
                return false;
            }

            return true;
        }

        VLIBEXP void replaceFactsInmemoryTable(std::string predicate,
                std::vector<std::vector<std::string>> &rows);

        // For JNI interface ...
        VLIBEXP void addInmemoryTable(std::string predicate,
                std::vector<std::vector<std::string>> &rows);

        VLIBEXP void addInmemoryTable(std::string predicate,
                PredId_t id, std::vector<std::vector<std::string>> &rows);

        //For RMFA check
        VLIBEXP void addInmemoryTable(PredId_t predicate,
                uint8_t arity,
                std::vector<uint64_t> &rows);

        VLIBEXP void addEDBTable(PredId_t predId,
                std::string tableType,
                std::shared_ptr<EDBTable> table);

        const EDBConf &getConf() const {
            return conf;
        }

        void setName(const std::string &name) {
            this->name = name;
        }

        const std::string &getName() const {
            return name;
        }

        //This method specifies whether the EDB source can change *during* reasoning
        //Most sources cannot, thus this method normally returns false
        bool canChange(PredId_t predId);

        //Some EDB sources do not allow variables (or constants) in some positions.
        //In this case, we cannot do a merge join but must apply sideway information
        //passing and other types of joins.
        bool isQueryAllowed(const Literal &query);

        //Some EDB sources require that the literal has no variable. For instance,
        //some string functions. In this case, we do not do any join
        bool acceptQueriesWithFreeVariables(const Literal &query);

        BuiltinFunction getBuiltinFunction(const Literal &query);

        void setContext(GBGraph *g, size_t step);

        void clearContext();

        ~EDBLayer() {
            for (int i = 0; i < tmpRelations.size(); ++i) {
                if (tmpRelations[i] != NULL) {
                    delete tmpRelations[i];
                }
            }
        }
};

#endif
