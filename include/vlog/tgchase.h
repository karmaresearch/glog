#ifndef _TG_CHASE_H
#define _TG_CHASE_H

#include <vlog/concepts.h>
#include <vlog/chase.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>
#include <vlog/tgsegment.h>

#include <vlog/gbgraph.h>
#include <vlog/gbruleexecutor.h>

#include <chrono>

struct CacheRetainEntry {
    size_t nnodes;
    std::shared_ptr<const TGSegment> seg;
};

class TGChase : public Chase {
    private:
        Program *program;
        std::vector<Rule> rules;
        EDBLayer &layer;

        GBGraph g;

        size_t currentIteration;
        PredId_t currentPredicate;
#ifdef WEBINTERFACE
        std::string currentRule;
#endif

        bool trackProvenance;
        std::vector<int> stratification;
        int nStratificationClasses;

        const bool cacheRetainEnabled;
        std::map<PredId_t, CacheRetainEntry> cacheRetain;
        std::chrono::duration<double, std::milli> durationRetain;

        GBRuleExecutor executor;

        bool executeRule(GBRuleInput &node);

        void createNewNodesWithProv(
                size_t ruleIdx, size_t step,
                std::shared_ptr<const TGSegment> seg,
                std::vector<std::shared_ptr<Column>> &provenance);

        std::shared_ptr<const TGSegment> retainVsNodeFast(
                std::shared_ptr<const TGSegment> existuples,
                std::shared_ptr<const TGSegment> newtuples);

        std::shared_ptr<const TGSegment> retain(
                PredId_t pred,
                std::shared_ptr<const TGSegment> newtuples);

    public:
        VLIBEXP TGChase(EDBLayer &layer, Program *program, bool useCacheRetain = true);

        VLIBEXP void run();

        Program *getProgram();

        EDBLayer &getEDBLayer();

        size_t getSizeTable(const PredId_t predid) const;

        FCIterator getTableItr(const PredId_t predid);

        FCTable *getTable(const PredId_t predid);

        size_t getCurrentIteration();

        size_t getNDerivedFacts();

        size_t getNnodes();

#ifdef WEBINTERFACE
        std::string getCurrentRule();

        PredId_t getCurrentPredicate();
#endif

};

#endif
