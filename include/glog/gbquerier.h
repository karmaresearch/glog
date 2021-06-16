#ifndef _GBQUERIER_H
#define _GBQUERIER_H

#include <trident/utils/json.h>
#include <glog/gbgraph.h>
#include <glog/gbsegment.h>
#include <vlog/concepts.h>
#include <vlog/edb.h>

#include <string>

class GBQuerier {
    private:
        const GBGraph &g;
        Program &p;
        EDBLayer &l;

        Literal getFact(PredId_t predId,
                std::shared_ptr<const TGSegment> data, size_t factId);

        std::string getTupleIDs(Literal &l);

        void exportNode(JSON &out, size_t nodeId, size_t factId);

        void exportNode(
                JSON &out,
                size_t nodeId, size_t factId,
                PredId_t nodePred,
                std::shared_ptr<const TGSegment> data,
                size_t ruleIdx,
                size_t step,
                const std::vector<size_t> &incomingEdges);

        void exportEDBNode(JSON &out, Literal &l, size_t factId);

    public:
        GBQuerier(const GBGraph &g, Program &p, EDBLayer &l) : g(g), p(p), l(l) {}

        JSON getDerivationTree(size_t nodeId, size_t factId);

        std::vector<std::string> getListPredicates() const;

        JSON getNodeDetailsWithPredicate(std::string predName) const;

        JSON getNodeFacts(size_t nodeId) const;

        bool checkSoundnessDerivationTree(JSON &root);
};

#endif
