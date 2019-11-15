#ifndef _TG_CHASE_H
#define _TG_CHASE_H

#include <vlog/concepts.h>
#include <vlog/chase.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>
#include <vlog/tgsegment.h>

#include <map>
#include <chrono>

struct TGChase_Node {
    size_t ruleIdx;
    size_t step;
    std::vector<size_t> incomingEdges;
    std::shared_ptr<const TGSegment> data;
    TGChase_Node() : ruleIdx(0), step(0) {}
};

struct TGChase_SuperNode {
    size_t ruleIdx;
    size_t step;
    std::vector<std::vector<size_t>> incomingEdges;
    TGChase_SuperNode() : ruleIdx(0), step(0) {}
};

class TGChase : public Chase {
    private:
        Program *program;
        std::vector<Rule> rules;
        EDBLayer &layer;
        size_t currentIteration;

        bool trackProvenance;

        std::map<PredId_t, std::shared_ptr<EDBTable>> edbTables;

#ifdef WEBINTERFACE
        PredId_t currentPredicate;
        std::string currentRule;
#endif

        std::map<PredId_t, std::vector<size_t>> pred2Nodes;
        std::vector<TGChase_Node> nodes;

        std::chrono::duration<double, std::milli> durationFirst;
        std::chrono::duration<double, std::milli> durationMergeSort;
        std::chrono::duration<double, std::milli> durationJoin;
        std::chrono::duration<double, std::milli> durationRetain;
        std::chrono::duration<double, std::milli> durationCreateHead;

        //Methods to execute the rule
        bool executeRule(TGChase_SuperNode &node);

        std::shared_ptr<const TGSegment> projectHead(const Literal &head,
                std::vector<size_t> &vars,
                std::shared_ptr<const TGSegment> intermediateResults,
                bool shouldSort,
                bool shouldDelDupl);

        void postprocessJoin(
            std::shared_ptr<const TGSegment> &intermediateResults,
            std::vector<std::shared_ptr<Column>> &intermediateResultsNodes,
            bool replace);

        int cmp(std::unique_ptr<TGSegmentItr> &inputLeft,
                std::unique_ptr<TGSegmentItr> &inputRight,
                std::pair<int, int> &joinVarPos);

        int cmp(std::unique_ptr<TGSegmentItr> &inputLeft,
                std::unique_ptr<TGSegmentItr> &inputRight);

        void computeVarPos(std::vector<size_t> &varsIntermediate,
                int bodyAtomIdx,
                const std::vector<Literal> &bodyAtoms,
                const Literal &head,
                std::pair<int, int> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight);

        std::shared_ptr<const TGSegment> processFirstAtom_EDB(
                const Literal &atom,
                std::vector<int> &copyVarPos);

        std::shared_ptr<const TGSegment> processFirstAtom_IDB(
                std::shared_ptr<const TGSegment> &input,
                std::vector<int> &copyVarPos,
                size_t nodeId);

        std::shared_ptr<const TGSegment> mergeNodes(
                std::vector<size_t> &nodeIdxs,
                std::vector<int> &copyVarPos);

        void mergejoin(
                std::shared_ptr<const TGSegment> inputLeft,
                std::shared_ptr<const TGSegment> inputRight,
                std::pair<int, int> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<SegmentInserter> &output);

        void join(
                std::shared_ptr<const TGSegment> inputLeft,
                std::vector<size_t> &bodyIdxs,
                std::pair<int, int> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<SegmentInserter> &output);

        void shouldSortDelDupls(const Literal &head,
                const std::vector<Literal> &bodyAtoms,
                const std::vector<std::vector<size_t>> &bodyNodes,
                bool &shouldSort,
                bool &shouldDelDupl);

        std::shared_ptr<const TGSegment> retain(
                PredId_t pred,
                std::shared_ptr<const TGSegment> newtuples);

        std::shared_ptr<const TGSegment> retainVsNodeFast(
                std::shared_ptr<const TGSegment> existuples,
                std::shared_ptr<const TGSegment> newtuples);

    public:
        VLIBEXP TGChase(EDBLayer &layer, Program *program);

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
