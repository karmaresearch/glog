#ifndef _GB_CHASE_H
#define _GB_CHASE_H

#include <vlog/concepts.h>
#include <vlog/chase.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>

#include <vlog/gbchase/gbsegment.h>
#include <vlog/gbchase/gbgraph.h>
#include <vlog/gbchase/gbruleexecutor.h>

#include <chrono>

typedef enum { GBCHASE, TGCHASE_STATIC, TGCHASE_DYNAMIC,TGCHASE_DYNAMIC_FULLPROV, PROBTGCHASE } GBChaseAlgorithm;

class GBChase : public Chase {
    protected:
        const GBGraph::ProvenanceType provenanceType;
        const bool filterQueryCont;
        const bool edbCheck;
        Program *program;
        EDBLayer &layer; //Stores the input data
        GBGraph g; //Stores the derivations
        std::vector<Rule> rules;

    private:
        std::vector<int> stratification;
        int nStratificationClasses;
        GBRuleExecutor executor; //Object that executes rules
        size_t currentIteration;

        PredId_t currentPredicate;
#ifdef WEBINTERFACE
        std::string currentRule;
#endif
        std::map<PredId_t, std::shared_ptr<FCTable>> cacheFCTables;
        std::set<PredId_t> predToBeRetainedEndStep;

        //Used for statistics
        size_t triggers;
        std::chrono::duration<double, std::milli> durationPreparation;
        std::chrono::duration<double, std::milli> durationRuleExec;
        //std::chrono::duration<double, std::milli> durationDebug;

        void prepareRuleExecutionPlans_queryContainment(
                std::vector<GBRuleInput> &newnodes,
                std::vector<std::vector<size_t>> &acceptableNodes,
                const size_t ruleIdx,
                const size_t step);

        void prepareRuleExecutionPlans_SNE(
                const std::vector<std::pair<bool,std::vector<size_t>>> nodesForRule,
                size_t ruleIdx,
                size_t step,
                size_t prevstep,
                std::vector<GBRuleInput> &newnodes);

        bool shouldRetainAtEnd(PredId_t pred);

    protected:
        bool shouldTrackProvenance() const {
            return provenanceType != GBGraph::ProvenanceType::NOPROV;
        }

        std::pair<bool, size_t> determineAdmissibleRule(
                const size_t &ruleIdx,
                const size_t stratumLevel,
                const size_t stepStratum,
                const size_t step) const;

        void determineAdmissibleRules(
                const std::vector<size_t> &ruleIdxs,
                const size_t stratumLevel,
                const size_t stepStratum,
                const size_t step,
                std::vector<std::pair<size_t,size_t>> &admissibleRules) const;

        void prepareRuleExecutionPlans(
                const size_t &ruleIdx,
                const size_t prevstep,
                const size_t step,
                std::vector<GBRuleInput> &newnodes);

        bool executeRule(GBRuleInput &node, bool cleanDuplicates = true);

        virtual size_t executeRulesInStratum(
                const std::vector<size_t> &ruleIdxs,
                const size_t stratumLevel,
                const size_t stepStratum,
                size_t &step);

    public:
        VLIBEXP GBChase(EDBLayer &layer, Program *program,
                bool useCacheRetain = true,
                GBGraph::ProvenanceType provenanceType = GBGraph::ProvenanceType::NOPROV,
                bool filterQueryCont = false,
                bool edbCheck = false,
                bool rewriteCliques = false);

        VLIBEXP virtual void run();

        Program *getProgram();

        EDBLayer &getEDBLayer();

        size_t getSizeTable(const PredId_t predid) const;

        FCIterator getTableItr(const PredId_t predid);

        FCTable *getTable(const PredId_t predid);

        size_t getCurrentIteration();

        size_t getNDerivedFacts();

        size_t getNnodes();

        size_t getNedges();

        size_t getNTriggers();

#ifdef WEBINTERFACE
        std::string getCurrentRule();

        PredId_t getCurrentPredicate();
#endif

};

#endif
