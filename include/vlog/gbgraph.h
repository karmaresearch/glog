#ifndef _GB_GRAPH_H
#define _GB_GRAPH_H

#include <vlog/concepts.h>
#include <vlog/tgsegment.h>

#include <map>

struct TGChase_Node {
    size_t ruleIdx;
    size_t step;
    std::vector<size_t> incomingEdges;
    std::shared_ptr<const TGSegment> data;
    TGChase_Node() : ruleIdx(0), step(0) {}
};

class GBGraph {
    private:
        std::map<PredId_t, std::vector<size_t>> pred2Nodes;
        std::vector<TGChase_Node> nodes;

    public:

        size_t getNNodes() const {
            return nodes.size();
        }

        const TGChase_Node getNode(size_t nodeId) const {
            return nodes[nodeId];
        }

        std::shared_ptr<const TGSegment> getNodeData(size_t nodeId) const {
            return nodes[nodeId].data;
        }

        bool areNodesWithPredicate(PredId_t predid) const {
            return pred2Nodes.count(predid);
        }

        const std::vector<TGChase_Node> &getNodes() const {
            return nodes;
        }

        const std::vector<size_t> &getNodeIDsWithPredicate(PredId_t predid) const {
            return pred2Nodes.at(predid);
        }

        void addNode(PredId_t predid, size_t ruleIdx, size_t step, std::shared_ptr<const TGSegment> data) {
            auto nodeId = getNNodes();
            nodes.emplace_back();
            TGChase_Node &outputNode = nodes.back();
            outputNode.ruleIdx = ruleIdx;
            outputNode.step = step;
            outputNode.data = data;
            pred2Nodes[predid].push_back(nodeId);
        }

        size_t getNDerivedFacts() const {
            size_t nderived = 0;
            for(auto &node : nodes) {
                nderived += node.data->getNRows();
            }
            return nderived;
        }

};

#endif
