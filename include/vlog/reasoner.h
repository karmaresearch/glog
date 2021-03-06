#ifndef REASONER_H
#define REASONER_H

#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/fctable.h>
#include <vlog/seminaiver.h>
#include <vlog/seminaiver_trigger.h>
#include <glog/gbchase.h>
#include <glog/tgchase.h>
#include <glog/probgbchase.h>
#include <vlog/consts.h>

#include <trident/kb/kb.h>
#include <trident/kb/querier.h>

#include <trident/sparql/query.h>
#include <vector>

#define QUERY_MAT 0
#define QUERY_ONDEM 1

typedef enum {TOPDOWN, MAGIC} ReasoningMode;

struct Metrics {
    size_t estimate;
    size_t intermediateResults;
    int countIntermediateQueries;
    double cost;
    int countRules;
    int countUniqueRules;
    int countIDBPredicates;
    uint8_t boundedness;
};

class Reasoner {
    private:

        const uint64_t threshold;

        void cleanBindings(std::vector<Term_t> &bindings,
                std::vector<uint8_t> * posJoins,
                TupleTable *input);

        FCBlock getBlockFromQuery(Literal constantsQuery, Literal &boundQuery,
                std::vector<uint8_t> *posJoins,
                std::vector<Term_t> *possibleValuesJoins);


        TupleIterator *getIncrReasoningIterator(Literal &query,
                std::vector<uint8_t> * posJoins,
                std::vector<Term_t> *possibleValuesJoins,
                EDBLayer &layer, Program &program,
                bool returnOnlyVars,
                std::vector<uint8_t> *sortByFields);

    public:

        Reasoner(const uint64_t threshold) : threshold(threshold) {}

        size_t estimate(Literal &query, std::vector<uint8_t> *posBindings,
                std::vector<Term_t> *valueBindings, EDBLayer &layer,
                Program &program);

        void getMetrics(Literal &query,
                std::vector<uint8_t> *posBindings,
                std::vector<Term_t> *valueBindings,
                EDBLayer &layer,
                Program &program,
                Metrics &metrics,
                int maxDepth,
                string& idbFeatures);

        VLIBEXP ReasoningMode chooseMostEfficientAlgo(Literal &query,
                EDBLayer &layer, Program &program,
                std::vector<uint8_t> *posBindings,
                std::vector<Term_t> *valueBindings);

        VLIBEXP TupleIterator *getIterator(Literal &query,
                std::vector<uint8_t> * posJoins,
                std::vector<Term_t> *possibleValuesJoins,
                EDBLayer &layer, Program &program,
                bool returnOnlyVars,
                std::vector<uint8_t> *sortByFields);

        VLIBEXP TupleIterator *getTopDownIterator(Literal &query,
                std::vector<uint8_t> * posJoins,
                std::vector<Term_t> *possibleValuesJoins,
                EDBLayer &layer, Program &program,
                bool returnOnlyVars,
                std::vector<uint8_t> *sortByFields);

        VLIBEXP TupleIterator *getMaterializationIterator(Literal &query,
                std::vector<uint8_t> * posJoins,
                std::vector<Term_t> *possibleValuesJoins,
                EDBLayer &layer, Program &program,
                bool returnOnlyVars,
                std::vector<uint8_t> *sortByFields);

        VLIBEXP TupleIterator *getIteratorWithMaterialization(SemiNaiver *sn,
                Literal &query,
                bool returnOnlyVars,
                std::vector<uint8_t> *sortByFields);

        VLIBEXP TupleIterator *getEDBIterator(Literal &query,
                std::vector<uint8_t> * posJoins,
                std::vector<Term_t> *possibleValuesJoins,
                EDBLayer &layer,
                bool returnOnlyVars,
                std::vector<uint8_t> *sortByFields);

        VLIBEXP TupleIterator *getMagicIterator(Literal &query,
                std::vector<uint8_t> * posJoins,
                std::vector<Term_t> *possibleValuesJoins,
                EDBLayer &layer, Program &program,
                bool returnOnlyVars,
                std::vector<uint8_t> *sortByFields);

        VLIBEXP TupleIterator *getProbMagicIterator(Literal &query,
                EDBLayer &layer, Program &program,
                bool returnOnlyVars,
                bool optDelProofsStaticAnalysis,
                std::string profilerPath="",
                std::string storematPath="");

        VLIBEXP TupleIterator *getTGMagicIterator(Literal &query,
                EDBLayer &layer, Program &program,
                bool returnOnlyVars,
                std::string profilerPath="",
                std::string storematPath="");

        VLIBEXP static std::shared_ptr<SemiNaiver> getSemiNaiver(EDBLayer &layer,
                Program *p, bool opt_intersect, bool opt_filtering, bool opt_threaded,
                TypeChase typeChase,
                int nthreads, int interRuleThreads, bool shuffleRules,
                Program *RMFC_check = NULL,
                std::string sameasAlgo = "");

        VLIBEXP static std::shared_ptr<TriggerSemiNaiver> getTriggeredSemiNaiver(
                EDBLayer &layer,
                Program *p,
                TypeChase typeChase);

        VLIBEXP static std::shared_ptr<GBChase> getGBChase(
                EDBLayer &layer,
                Program *p,
                GBChaseAlgorithm typeChase = GBChaseAlgorithm::GBCHASE,
                bool queryCont = true,
                bool edbCheck = true,
                bool rewriteCliques = true,
                std::string param1 = "");

        VLIBEXP static std::shared_ptr<GBChase> getProbTGChase(
                EDBLayer &layer,
                Program *p,
                bool optDelProofsStaticAnalysis);

        int getNumberOfIDBPredicates(Literal&, Program&);

        ~Reasoner() {
        }
};
#endif
