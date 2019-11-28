#ifndef _SEMINAIVER_H
#define _SEMINAIVER_H

#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/fctable.h>
#include <vlog/ruleexecplan.h>
#include <vlog/ruleexecdetails.h>
#include <vlog/chasemgmt.h>
#include <vlog/consts.h>
#include <vlog/chase.h>

#include <trident/model/table.h>

#include <vector>
#include <unordered_map>

struct StatIteration {
    size_t iteration;
    const Rule *rule;
    double time;
    bool derived;

    bool operator <(const StatIteration &it) const {
        return time > it.time;
    }
};

typedef std::unordered_map<std::string, FCTable*> EDBCache;
class ResultJoinProcessor;
class SemiNaiver: public Chase {
    private:
        std::vector<RuleExecutionDetails> allEDBRules;
        bool opt_intersect;
        bool opt_filtering;
        bool multithreaded;

        std::vector<FCBlock> listDerivations;

        bool ignoreDuplicatesElimination;
        std::vector<int> stratification;
        int nStratificationClasses;
        Program *RMFC_program;

#ifdef WEBINTERFACE
        std::string currentRule;
        PredId_t currentPredicate;
#endif

    private:
        FCIterator getTableFromIDBLayer(const Literal & literal,
                const size_t minIteration,
                TableFilterer *filter);

        FCIterator getTableFromIDBLayer(const Literal & literal,
                const size_t minIteration,
                const size_t maxIteration,
                TableFilterer *filter);

        size_t countAllIDBs();

        bool bodyChangedSince(Rule &rule, size_t iteration);

        bool checkIfAtomsAreEmpty(const RuleExecutionDetails &ruleDetails,
                const RuleExecutionPlan &plan,
                size_t limitView,
                std::vector<size_t> &cards);

        void processRuleFirstAtom(const uint8_t nBodyLiterals,
                const Literal *bodyLiteral,
                std::vector<Literal> &heads,
                const size_t min,
                const size_t max,
                int &processedTables,
                const bool lastLiteral,
                const size_t iteration,
                const RuleExecutionDetails &ruleDetails,
                const uint8_t orderExecution,
                std::vector<std::pair<uint8_t, uint8_t>> *filterValueVars,
                ResultJoinProcessor *joinOutput);

        void reorderPlan(RuleExecutionPlan &plan,
                const std::vector<size_t> &cards,
                const std::vector<Literal> &headLiteral,
                bool copyAllVars);

        void reorderPlanForNegatedLiterals(RuleExecutionPlan &plan,
                const std::vector<Literal> &heads);

        bool executeRules(std::vector<RuleExecutionDetails> &allEDBRules,
                std::vector<std::vector<RuleExecutionDetails>> &allIDBRules,    // one entry for each stratification class
                std::vector<StatIteration> &costRules,
                const size_t limitView,
                bool fixpoint, unsigned long *timeout = NULL);

        bool executeRule(RuleExecutionDetails &ruleDetails,
                std::vector<Literal> &heads,
                const size_t iteration,
                const size_t limitView,
                std::vector<ResultJoinProcessor*> *finalResultContainer);

        size_t estimateCardTable(const Literal &literal,
                const size_t minIteration,
                const size_t maxIteration);

    protected:
        TypeChase typeChase;
        bool checkCyclicTerms;
        bool foundCyclicTerms;
        PredId_t predIgnoreBlock; //RMSA
        bool ignoreExistentialRules;
        std::shared_ptr<ChaseMgmt> chaseMgmt;

        std::vector<FCTable *>predicatesTables;
        EDBLayer &layer;
        Program *program;
        std::vector<std::vector<RuleExecutionDetails>> allIDBRules; // one entry for each stratification class
        size_t iteration;
        int nthreads;
        uint64_t triggers;

        bool executeRule(RuleExecutionDetails &ruleDetails,
                const size_t iteration,
                const size_t limitView,
                std::vector<ResultJoinProcessor*> *finalResultContainer);

        virtual FCIterator getTableFromEDBLayer(const Literal & literal);

        virtual long getNLastDerivationsFromList();

        virtual void saveDerivationIntoDerivationList(FCTable *endTable);

        virtual bool executeUntilSaturation(
                std::vector<RuleExecutionDetails> &ruleset,
                std::vector<StatIteration> &costRules,
                size_t limitView,
                bool fixpoint, unsigned long *timeout = NULL);

        void prepare(size_t lastExecution, int singleRuleToCheck, std::vector<RuleExecutionDetails> &allrules);

        void setIgnoreDuplicatesElimination() {
            ignoreDuplicatesElimination = true;
        }

    public:
        VLIBEXP SemiNaiver(EDBLayer &layer,
                Program *program, bool opt_intersect,
                bool opt_filtering, bool multithreaded,
                TypeChase chase, int nthreads, bool shuffleRules,
                bool ignoreExistentialRule, Program *RMFC_check = NULL);

        //disable restricted chase
        VLIBEXP SemiNaiver(EDBLayer &layer,
                Program *program, bool opt_intersect,
                bool opt_filtering, bool multithreaded,
                int nthreads, bool shuffleRules,
                bool ignoreExistentialRules) :
            SemiNaiver(layer, program, opt_intersect, opt_filtering,
                    multithreaded, TypeChase::SKOLEM_CHASE, nthreads, shuffleRules,
                    ignoreExistentialRules) {
            }

        Program *get_RMFC_program() {
            return RMFC_program;
        }

        bool opt_filter() {
            return opt_filtering;
        }

        bool opt_inter() {
            return opt_intersect;
        }

        std::shared_ptr<ChaseMgmt> getChaseManager() {
            return chaseMgmt;
        }

        virtual FCTable *getTable(const PredId_t pred, const uint8_t card);

        VLIBEXP FCTable *getTable(const PredId_t predid);

        VLIBEXP void run() {
            run(NULL, false);
        }

        VLIBEXP void run(unsigned long *timeout,
                bool checkCyclicTerms = false) {
            run(0, 1, timeout, checkCyclicTerms, -1, -1);
        }

        VLIBEXP void run(size_t lastIteration,
                size_t iteration,
                unsigned long *timeout = NULL,
                bool checkCyclicTerms = false,
                int singleRule = -1,
                PredId_t predIgnoreBlock = -1);

        std::ostream& dumpTables(std::ostream &os) {
            for (PredId_t i = 0; i < MAX_NPREDS; ++i) {
                FCTable *table = predicatesTables[i];
                if (table != NULL && !table->isEmpty()) {
                    char buffer[MAX_TERM_SIZE];

                    os << "Table " << getProgram()->getPredicateName(i) << std::endl;
                    FCIterator itr = table->read(0);
                    const uint8_t sizeRow = table->getSizeRow();
                    while (!itr.isEmpty()) {
                        std::shared_ptr<const FCInternalTable> t = itr.getCurrentTable();
                        FCInternalTableItr *iitr = t->getIterator();
                        while (iitr->hasNext()) {
                            iitr->next();
                            std::string row = "    ";
                            row += std::to_string(iitr->getCurrentIteration());
                            for (uint8_t m = 0; m < sizeRow; ++m) {
                                row += "\t" + std::to_string(iitr->getCurrentValue(m));
                            }
                            os << row << std::endl;
                        }
                        t->releaseIterator(iitr);
                        itr.moveNextCount();
                    }
                }
            }

            return os;
        }

        FCIterator getTable(const Literal &literal, const size_t minIteration,
                const size_t maxIteration) {
            return getTable(literal, minIteration, maxIteration, NULL);
        }

        VLIBEXP FCIterator getTableItr(const PredId_t predid);

        size_t getSizeTable(const PredId_t predid) const;

        bool isEmpty(const PredId_t predid) const;

        std::vector<FCBlock> &getDerivationsSoFar() {
            return listDerivations;
        }

        VLIBEXP void createGraphRuleDependency(std::vector<int> &nodes,
                std::vector<std::pair<int, int>> &edges);

        Program *getProgram() {
            return program;
        }

        EDBLayer &getEDBLayer() {
            return layer;
        }

        bool isFoundCyclicTerms() {
            return foundCyclicTerms;
        }

        void addDataToIDBRelation(const Predicate pred, FCBlock block);

        size_t estimateCardinality(const Literal &literal, const size_t min,
                const size_t max);

        virtual ~SemiNaiver();

        static std::pair<uint8_t, uint8_t> removePosConstants(
                std::pair<uint8_t, uint8_t> columns,
                const Literal &literal);

        virtual FCIterator getTable(const Literal &literal, const size_t minIteration,
                const size_t maxIteration, TableFilterer *filter);

        void checkAcyclicity(int singleRule = -1, PredId_t predIgnoreBlock = -1) {
            run(0, 1, NULL, true, singleRule, predIgnoreBlock);
        }

        //Statistics methods
        VLIBEXP void printCountAllIDBs(std::string prefix);

        size_t getCurrentIteration();

#ifdef WEBINTERFACE
        std::string getCurrentRule();

        PredId_t getCurrentPredicate();
#endif

};

#endif
