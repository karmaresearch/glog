#ifndef _GB_RULEEXECUTOR_H
#define _GB_RULEEXECUTOR_H

#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>

#include <glog/gbgraph.h>
#include <glog/gbsegment.h>
#include <glog/gbsegmentinserter.h>

#include <chrono>

struct GBRuleInput {
    size_t ruleIdx;
    size_t step;
    bool retainFree;
    std::vector<std::vector<size_t>> incomingEdges;
    GBRuleInput() : ruleIdx(0), step(0), retainFree(false) {}
};

struct GBRuleOutput {
    std::shared_ptr<const TGSegment> segment;
    std::vector<std::shared_ptr<Column>> nodes;
    bool uniqueTuples;
    GBRuleOutput() : uniqueTuples(false) {}
};

typedef enum {DUR_FIRST, DUR_MERGE, DUR_JOIN, DUR_HEAD //,DUR_PREP2TO1
} DurationType;
typedef enum { N_BDY_ATOMS } StatType;


#define N_ATTEMPTS_ENABLE_DUPL_DEL 5
class DuplManager {
    private:
        bool enabled;
        GBSegmentInserterEntities l;
        GBSegmentInserterEntities r;

    public:
        DuplManager(GBSegmentInserter *output);
        bool left(const Term_t &t) const;
        bool right(const Term_t &t) const;
};

class GBRuleExecutor {
    private:
        std::chrono::duration<double, std::milli> durationFirst;
        std::chrono::duration<double, std::milli> durationMergeSort;
        std::chrono::duration<double, std::milli> durationJoin;
        std::chrono::duration<double, std::milli> durationCreateHead;

        std::chrono::duration<double, std::milli> lastDurationFirst;
        std::chrono::duration<double, std::milli> lastDurationMergeSort;
        std::chrono::duration<double, std::milli> lastDurationJoin;
        std::chrono::duration<double, std::milli> lastDurationCreateHead;
        std::string bdyAtoms;

        Program *program; //used only for debugging purposes
        EDBLayer &layer;
        std::map<PredId_t, std::shared_ptr<EDBTable>> edbTables;

        const GBGraph::ProvenanceType provenanceType;
        GBGraph &g;
        std::vector<size_t> noBodyNodes;

        void shouldSortAndRetainEDBSegments(
                bool &shouldSort,
                bool &shouldRetainUnique,
                const Literal &atom,
                std::vector<int> &copyVarPos);

        std::shared_ptr<const TGSegment> projectTuples(
                std::shared_ptr<const TGSegment> tuples,
                const std::vector<int> &posVariables);

        std::shared_ptr<const TGSegment> projectTuples_structuresharing(
                std::shared_ptr<const TGSegment> tuples,
                const std::vector<int> &posVariables,
                bool isSorted);

        std::shared_ptr<const TGSegment> projectHead(const Literal &head,
                std::vector<size_t> &vars,
                std::shared_ptr<const TGSegment> intermediateResults,
                bool shouldSort,
                bool shouldDelDupl,
                std::vector<std::shared_ptr<Column>> &intermediateResultsNodes,
                GBRuleInput &node);

        std::shared_ptr<const TGSegment> addConstants(
                std::shared_ptr<const TGSegment> tg,
                std::vector<std::pair<size_t, Term_t>> constants);

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

        std::shared_ptr<const TGSegment> processAtom_EDB(
                const Literal &atom,
                std::vector<int> &copyVarPos);

        std::shared_ptr<const TGSegment> processAtom_IDB(
                const Literal &atom,
                std::vector<size_t> &nodeIdxs,
                std::vector<int> &copyVarPos,
                bool lazyMode,
                bool replaceOffsets);

        void joinTwoOne_EDB(
                std::shared_ptr<const TGSegment> inputLeft,
                std::shared_ptr<const TGSegment> inputRight,
                int joinLeftVarPos,
                std::vector<int> &copyVarPosLeft,
                std::unique_ptr<GBSegmentInserter> &output);

        void joinTwoOne(
                std::shared_ptr<const TGSegment> inputLeft,
                std::shared_ptr<const TGSegment> inputRight,
                int joinLeftVarPos,
                std::vector<int> &copyVarPosLeft,
                const int copyNodes,
                std::unique_ptr<GBSegmentInserter> &output);

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

        bool shouldTrackProvenance() {
            return provenanceType != GBGraph::ProvenanceType::NOPROV;
        }

        SegProvenanceType getSegProvenanceType() const;

    public:
        GBRuleExecutor(GBGraph &g, EDBLayer &layer,
                Program *program) :
            durationMergeSort(0),
            durationJoin(0),
            durationCreateHead(0),
            durationFirst(0),
            lastDurationMergeSort(0),
            lastDurationJoin(0),
            lastDurationCreateHead(0),
            lastDurationFirst(0),
            bdyAtoms(""),
            program(program),
            provenanceType(g.getProvenanceType()),
            g(g), layer(layer)
    {
    }

        std::vector<GBRuleOutput> executeRule(Rule &rule, GBRuleInput &node);

        std::chrono::duration<double, std::milli> getDuration(DurationType typ);

        std::string getStat(StatType typ);

        void printStats();
};

#endif
