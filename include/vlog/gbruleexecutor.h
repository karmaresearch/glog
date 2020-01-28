#ifndef _GB_RULEEXECUTOR_H
#define _GB_RULEEXECUTOR_H

#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>

#include <vlog/gbgraph.h>
#include <vlog/tgsegment.h>

#include <chrono>

struct GBRuleInput {
    size_t ruleIdx;
    size_t step;
    std::vector<std::vector<size_t>> incomingEdges;
    GBRuleInput() : ruleIdx(0), step(0) {}
};

typedef std::pair<std::shared_ptr<const TGSegment>,
        std::vector<std::shared_ptr<Column>>> OutputRule;

class GBRuleExecutor {
    private:
        std::chrono::duration<double, std::milli> durationFirst;
        std::chrono::duration<double, std::milli> durationMergeSort;
        std::chrono::duration<double, std::milli> durationJoin;
        std::chrono::duration<double, std::milli> durationCreateHead;

        EDBLayer &layer;
        std::map<PredId_t, std::shared_ptr<EDBTable>> edbTables;

        const bool trackProvenance;
        GBGraph &g;
        std::vector<size_t> noBodyNodes;

        std::shared_ptr<const TGSegment> projectHead(const Literal &head,
                std::vector<size_t> &vars,
                std::shared_ptr<const TGSegment> intermediateResults,
                bool shouldSort,
                bool shouldDelDupl);

        std::shared_ptr<const Segment> postprocessJoin(
                std::shared_ptr<const Segment> &intermediateResults,
                std::vector<std::shared_ptr<Column>> &intermediateResultsNodes);

        void computeVarPos(std::vector<size_t> &varsIntermediate,
                int bodyAtomIdx,
                const std::vector<Literal> &bodyAtoms,
                const Literal &head,
                std::vector<std::pair<int, int>> &joinVarPos,
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
                const std::vector<size_t> &nodesLeft,
                std::shared_ptr<const TGSegment> inputRight,
                const std::vector<size_t> &nodesRight,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<SegmentInserter> &output);

        void leftjoin(
                std::shared_ptr<const TGSegment> inputLeft,
                std::vector<size_t> &bodyNodeIdxs,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::unique_ptr<SegmentInserter> &output);

        void join(
                std::shared_ptr<const TGSegment> inputLeft,
                const std::vector<size_t> &nodesLeft,
                std::vector<size_t> &bodyIdxs,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<SegmentInserter> &output);

        void shouldSortDelDupls(const Literal &head,
                const std::vector<Literal> &bodyAtoms,
                const std::vector<std::vector<size_t>> &bodyNodes,
                bool &shouldSort,
                bool &shouldDelDupl);

    public:
        GBRuleExecutor(bool trackProvenance, GBGraph &g, EDBLayer &layer) :
            durationMergeSort(0),
            durationJoin(0),
            durationCreateHead(0),
            durationFirst(0),
            trackProvenance(trackProvenance),
            g(g), layer(layer) {
            }

        static std::shared_ptr<const TGSegment> fromSeg2TGSeg(
                std::shared_ptr<const Segment> seg,
                size_t nodeId, bool isSorted, uint8_t sortedField,
                bool trackProvenance);

        OutputRule executeRule(Rule &rule, GBRuleInput &node);

        void printStats();
};

#endif
