#ifndef _GB_CHASE_H
#define _GB_CHASE_H

#include <vlog/concepts.h>
#include <vlog/chase.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>
#include <vlog/tgsegment.h>

#include <vlog/gbgraph.h>
#include <vlog/gbruleexecutor.h>

#include <chrono>

class GBChase : public Chase {
    private:
        bool trackProvenance;
        Program *program;
        std::vector<Rule> rules;
        std::vector<int> stratification;
        int nStratificationClasses;

        EDBLayer &layer; //Stores the input data
        GBGraph g; //Stores the derivations
        GBRuleExecutor executor; //Object that executes rules
        size_t currentIteration;

        PredId_t currentPredicate;
#ifdef WEBINTERFACE
        std::string currentRule;
#endif

        bool executeRule(GBRuleInput &node);

        void createNewNodesWithProv(
                size_t ruleIdx, size_t step,
                std::shared_ptr<const TGSegment> seg,
                std::vector<std::shared_ptr<Column>> &provenance);

    public:
        VLIBEXP GBChase(EDBLayer &layer, Program *program, bool useCacheRetain = true);

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
