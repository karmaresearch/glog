#ifndef _GB_RULEEXECUTOR_H
#define _GB_RULEEXECUTOR_H

#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>

#include <vlog/gbchase/gbgraph.h>
#include <vlog/gbchase/gbsegment.h>
#include <vlog/gbchase/gbsegmentinserter.h>

#include <chrono>

struct GBRuleInput {
    size_t ruleIdx;
    size_t step;
    std::vector<std::vector<size_t>> incomingEdges;
    GBRuleInput() : ruleIdx(0), step(0) {}
};

struct GBRuleOutput {
    std::shared_ptr<const TGSegment> segment;
    std::vector<std::shared_ptr<Column>> nodes;
    bool uniqueTuples;

    GBRuleOutput() : uniqueTuples(false) {}
};


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

        std::shared_ptr<const TGSegment> projectTuples(
                std::shared_ptr<const TGSegment> tuples,
                const std::vector<int> &posKnownVariables);

        std::shared_ptr<const TGSegment> projectHead(const Literal &head,
                std::vector<size_t> &vars,
                std::shared_ptr<const TGSegment> intermediateResults,
                bool shouldSort,
                bool shouldDelDupl);

        void computeVarPos(std::vector<size_t> &varsIntermediate,
                int bodyAtomIdx,
                const std::vector<Literal> &bodyAtoms,
                const std::vector<Literal> &heads,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight);

        void addBuiltinFunctions(std::vector<BuiltinFunction> &out,
                const std::vector<Literal> &atoms,
                const std::vector<size_t> &vars);

        std::shared_ptr<const TGSegment> processFirstAtom_EDB(
                const Literal &atom,
                std::vector<int> &copyVarPos);

        void mergejoin(
                std::shared_ptr<const TGSegment> inputLeft,
                const std::vector<size_t> &nodesLeft,
                std::shared_ptr<const TGSegment> inputRight,
                const std::vector<size_t> &nodesRight,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<GBSegmentInserter> &output);

        void nestedloopjoin(
                std::shared_ptr<const TGSegment> inputLeft,
                const std::vector<size_t> &nodesLeft,
                std::shared_ptr<const TGSegment> inputRight,
                const std::vector<size_t> &nodesRight,
                const Literal &literalRight,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<GBSegmentInserter> &output);

        void leftjoin(
                std::shared_ptr<const TGSegment> inputLeft,
                const std::vector<size_t> &nodesLeft,
                std::shared_ptr<const TGSegment> inputRight,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::unique_ptr<GBSegmentInserter> &output,
                const bool copyOnlyLeftNode = false);

        void join(
                std::shared_ptr<const TGSegment> inputLeft,
                const std::vector<size_t> &nodesLeft,
                std::vector<size_t> &nodesRight,
                const Literal &literalRight,
                std::vector<std::pair<int, int>> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<GBSegmentInserter> &output);

        void shouldSortDelDupls(const Literal &head,
                const std::vector<Literal> &bodyAtoms,
                const std::vector<std::vector<size_t>> &bodyNodes,
                bool &shouldSort,
                bool &shouldDelDupl);

        std::shared_ptr<const TGSegment> performRestrictedCheck(Rule &rule,
                std::shared_ptr<const TGSegment> tuples,
                const std::vector<size_t> &varTuples);

        std::shared_ptr<const TGSegment> addExistentialVariables(
                Rule &rule,
                std::shared_ptr<const TGSegment> tuples,
                std::vector<size_t> &vars);

    public:
        GBRuleExecutor(bool trackProvenance, GBGraph &g, EDBLayer &layer) :
            durationMergeSort(0),
            durationJoin(0),
            durationCreateHead(0),
            durationFirst(0),
            trackProvenance(trackProvenance),
            g(g), layer(layer) {
            }

        /*static std::shared_ptr<const TGSegment> fromSeg2TGSeg(
                std::shared_ptr<const Segment> seg,
                size_t nodeId, bool isSorted, uint8_t sortedField,
                bool trackProvenance);*/

        std::vector<GBRuleOutput> executeRule(Rule &rule, GBRuleInput &node);

        void printStats();
};

#endif
