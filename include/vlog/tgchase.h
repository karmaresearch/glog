#ifndef _TG_CHASE_H
#define _TG_CHASE_H

#include <vlog/concepts.h>
#include <vlog/chase.h>
#include <vlog/edb.h>
#include <vlog/fcinttable.h>

#include <map>
#include <chrono>

struct TGChase_Node {
    size_t ruleIdx;
    size_t step;
    std::vector<size_t> incomingEdges;
    std::shared_ptr<const Segment> data;
    TGChase_Node() : ruleIdx(0), step(0) {}
};

class TGChase : public Chase {
    private:
        Program *program;
        std::vector<Rule> rules;
        EDBLayer &layer;
        size_t currentIteration;

        std::map<PredId_t, std::unique_ptr<FCInternalTable>> edbTables;

#ifdef WEBINTERFACE
        PredId_t currentPredicate;
        std::string currentRule;
#endif

        std::map<PredId_t, std::vector<size_t>> pred2Nodes;
        std::vector<TGChase_Node> nodes;

        std::chrono::duration<double, std::milli> durationFirst;
        std::chrono::duration<double, std::milli> durationJoin;
        std::chrono::duration<double, std::milli> durationRetain;
        std::chrono::duration<double, std::milli> durationCreateHead;

        //Methods to execute the rule
        bool executeRule(size_t nodeId);

        std::shared_ptr<const Segment> projectHead(const Literal &head,
                std::vector<size_t> &vars,
                std::shared_ptr<const Segment> intermediateResults);

        int cmp(std::unique_ptr<SegmentIterator> &inputLeft,
                std::unique_ptr<SegmentIterator> &inputRight,
                std::pair<int, int> &joinVarPos);

        int cmp(std::unique_ptr<SegmentIterator> &inputLeft,
                std::unique_ptr<SegmentIterator> &inputRight);

        void computeVarPos(std::vector<size_t> &varsIntermediate,
                int bodyAtomIdx,
                const std::vector<Literal> &bodyAtoms,
                const Literal &head,
                std::pair<int, int> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight);

        void processFirstAtom_EDB(
                const Literal &atom,
                std::vector<int> &copyVarPos,
                std::unique_ptr<SegmentInserter> &output);

        void processFirstAtom_IDB(
                std::shared_ptr<const Segment> &input,
                std::vector<int> &copyVarPos,
                std::unique_ptr<SegmentInserter> &output);

        void recursiveCreateNode(
                const size_t step,
                const size_t ruleIdx,
                std::vector<std::vector<size_t>> &input,
                std::vector<size_t> &currentRow,
                const size_t columnIdx);

        void mergejoin(
                std::shared_ptr<const Segment> inputLeft,
                std::shared_ptr<const Segment> inputRight,
                std::pair<int, int> &joinVarPos,
                std::vector<int> &copyVarPosLeft,
                std::vector<int> &copyVarPosRight,
                std::unique_ptr<SegmentInserter> &output);

        std::shared_ptr<const Segment> retain(
                PredId_t pred,
                std::shared_ptr<const Segment> newtuples);

        std::shared_ptr<const Segment> retainVsNode(
                std::shared_ptr<const Segment> existuples,
                std::shared_ptr<const Segment> newtuples);

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
