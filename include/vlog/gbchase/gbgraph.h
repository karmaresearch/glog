#ifndef _GB_GRAPH_H
#define _GB_GRAPH_H

#include <vlog/concepts.h>
#include <vlog/chasemgmt.h>

#include <vlog/gbchase/gbsegment.h>

#include <map>

class GBGraph {
    private:
        struct GBGraph_Node {
            size_t ruleIdx;
            size_t step;
            std::vector<size_t> incomingEdges;
            std::shared_ptr<const TGSegment> data;
            GBGraph_Node() : ruleIdx(0), step(0) {}
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

        std::shared_ptr<const TGSegment> getNodeData(size_t nodeId) const {
            return nodes[nodeId].data;
        }

        bool areNodesWithPredicate(PredId_t predid) const {
            return pred2Nodes.count(predid);
        }

        const std::vector<size_t> &getNodeIDsWithPredicate(PredId_t predid) const {
            return pred2Nodes.at(predid);
        }

        void addNode(PredId_t predid, size_t ruleIdx,
                size_t step, std::shared_ptr<const TGSegment> data);

        void replaceEqualTerms(std::shared_ptr<const TGSegment> data);

        size_t getNDerivedFacts() const {
            size_t nderived = 0;
            for(auto &node : nodes) {
                nderived += node.data->getNRows();
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

        void printStats() {
            LOG(INFOL) << "Time retain (ms): " << durationRetain.count();
        }
};

#endif
