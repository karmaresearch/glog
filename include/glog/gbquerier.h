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

        void getLeaves(
                size_t nodeId, size_t factId,
                PredId_t nodePred,
                std::shared_ptr<const TGSegment> data,
                size_t ruleIdx,
                size_t step,
                const std::vector<size_t> &incomingEdges,
                std::vector<Literal> &out);

        void exportEDBNode(JSON &out, Literal &l, size_t factId);

        void exportEDBNode(Literal &l, size_t factId, std::vector<Literal> &out);

        void getMappings(const Literal &l,
                const std::vector<uint64_t> &row,
                std::vector<std::pair<Term_t, Term_t>> &mappings);

        Literal ground(const Literal &l, const std::vector<
                std::pair<Term_t, Term_t>> &mappings,
                bool &ok);

    public:
        GBQuerier(const GBGraph &g, Program &p, EDBLayer &l) : g(g), p(p), l(l) {}

        JSON getDerivationTree(size_t nodeId, size_t factId);

        JSON getDerivationTree(
                std::shared_ptr<const TGSegment> data,
                size_t nodeId,
                size_t factId,
                PredId_t predId,
                size_t ruleIdx,
                size_t step,
                const std::vector<size_t> &incomingEdges);

        std::vector<Term_t> getLeavesInDerivationTree(
                size_t nodeId,
                size_t factId,
                std::vector<Literal> &out);

        std::vector<std::string> getListPredicates() const;

        std::string getTermText(Term_t t) const;

        JSON getNodeDetailsWithPredicate(std::string predName) const;

        JSON getNodeFacts(size_t nodeId) const;

        std::map<std::string, std::vector<std::vector<std::string>>> getAllFacts() const;

        bool checkSoundnessDerivationTree(JSON &root);

};

#endif
