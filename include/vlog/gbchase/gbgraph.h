#ifndef _GB_GRAPH_H
#define _GB_GRAPH_H

#include <vlog/concepts.h>
#include <vlog/chasemgmt.h>

#include <vlog/gbchase/gbsegment.h>

#include <map>

class GBGraph {
    private:
        struct GBGraph_Node {
            private:
                std::shared_ptr<const TGSegment> data;
            public:
                PredId_t predid;
                size_t ruleIdx;
                size_t step;
                std::vector<size_t> incomingEdges;

                GBGraph_Node() : ruleIdx(0), step(0) {}

                std::shared_ptr<const TGSegment> getData() const {
                    return data;
                }

                void setData(std::shared_ptr<const TGSegment> data) {
                    this->data = data;
                }
        };

        struct CacheRetainEntry {
            size_t nnodes;
            std::shared_ptr<const TGSegment> seg;
        };

        std::map<PredId_t, std::vector<size_t>> pred2Nodes;
        std::vector<GBGraph_Node> nodes;
        const bool trackProvenance;
        const bool cacheRetainEnabled;
        std::map<PredId_t, CacheRetainEntry> cacheRetain;
        std::chrono::duration<double, std::milli> durationRetain;
        uint64_t counterNullValues;

        std::shared_ptr<const TGSegment> retainVsNodeFast(
                std::shared_ptr<const TGSegment> existuples,
                std::shared_ptr<const TGSegment> newtuples);

    public:
        GBGraph(bool trackProvenance, bool cacheRetainEnabled) :
            trackProvenance(trackProvenance),
            cacheRetainEnabled(cacheRetainEnabled),
            durationRetain(0) {
                counterNullValues = RULE_SHIFT(1);
            }

        size_t getNNodes() const {
            return nodes.size();
        }

        size_t getNodeStep(size_t nodeId) const {
            return nodes[nodeId].step;
        }

        size_t getNodeSize(size_t nodeId) const {
            return nodes[nodeId].getData()->getNRows();
        }

        std::shared_ptr<const TGSegment> getNodeData(size_t nodeId) const {
            return nodes[nodeId].getData();
        }

        PredId_t getNodePredicate(size_t nodeId) const {
            return nodes[nodeId].predid;
        }

        bool areNodesWithPredicate(PredId_t predId) const {
            return pred2Nodes.count(predId);
        }

        const std::vector<size_t> &getNodeIDsWithPredicate(
                PredId_t predId) const {
            return pred2Nodes.at(predId);
        }

        std::shared_ptr<const TGSegment> mergeNodes(
                const std::vector<size_t> &nodeIdxs,
                std::vector<int> &copyVarPos) const;

        void addNode(PredId_t predId, size_t ruleIdx,
                size_t step, std::shared_ptr<const TGSegment> data);

        void replaceEqualTerms(
                size_t ruleIdx,
                size_t step,
                std::shared_ptr<const TGSegment> data);

        size_t getNDerivedFacts() const {
            size_t nderived = 0;
            for(auto &node : nodes) {
                nderived += node.getData()->getNRows();
            }
            return nderived;
        }

        uint64_t getCounterNullValues() const {
            return counterNullValues;
        }

        void setCounterNullValues(uint64_t c) {
            counterNullValues = c;
        }

        std::shared_ptr<const TGSegment> retain(
                PredId_t pred,
                std::shared_ptr<const TGSegment> newtuples);

        //Returns the number of retained tuples. The new node will get the last
        //step and will be assigned to rule ~0ul
        uint64_t mergeNodesWithPredicateIntoOne(PredId_t predId);

        void printStats() {
            LOG(INFOL) << "Time retain (ms): " << durationRetain.count();
        }
};

#endif
