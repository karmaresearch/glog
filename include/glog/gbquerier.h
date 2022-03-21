#ifndef _GBQUERIER_H
#define _GBQUERIER_H

#include <trident/utils/json.h>
#include <glog/gbgraph.h>
#include <glog/gbsegment.h>
#include <vlog/concepts.h>
#include <vlog/edb.h>

#include <string>

#define MAX_LENGTH 200
#define MAX_TUPLE_ARITY 3
class DuplicateChecker
{
    private:
        bool enabled;
        std::vector<std::pair<PredId_t, int>> predicates;
        std::vector<Term_t> tuples;
        size_t idxLastElement;

    public:
        DuplicateChecker() :
            enabled(false), predicates(MAX_LENGTH),
            tuples(MAX_TUPLE_ARITY * MAX_LENGTH),
            idxLastElement(0) {}
        bool isEnabled() const { return enabled;}
        void setEnabled() { enabled = true; }
        bool isRedundant(const Literal &l) const;
        bool isRedundant(const PredId_t &p,
                const std::vector<Term_t> &row) const;
        void pushFact(const Literal &l);
        void pushFact(const PredId_t &p, int arity,
                const std::vector<Term_t> &row);
        size_t getMark();
        void setMark(size_t mark);
};

class IncompleteProofInfo
{
    public:
        size_t start;
        std::vector<Literal> leaves;
        std::vector<std::pair<Term_t, Term_t>> mappings;
        size_t j;
};

class GBQuerier {
    private:
        const GBGraph &g;
        Program &p;
        EDBLayer &l;


        Literal getFact(PredId_t predId,
                std::shared_ptr<const TGSegment> data, size_t factId);

        std::string getTupleIDs(Literal &l);

        bool exportNode(JSON &out, size_t nodeId, size_t factId,
                DuplicateChecker *checker);

        bool exportNode(
                JSON &out,
                size_t nodeId, size_t factId,
                PredId_t nodePred,
                std::shared_ptr<const TGSegment> data,
                size_t ruleIdx,
                size_t step,
                const std::vector<size_t> &incomingEdges,
                DuplicateChecker *checker);

        bool getLeaves(
                size_t nodeId, size_t factId,
                PredId_t nodePred,
                std::shared_ptr<const TGSegment> data,
                size_t ruleIdx,
                size_t step,
                const std::vector<size_t> &incomingEdges,
                std::vector<std::vector<Literal>> &out,
                DuplicateChecker *checker);

        //This method is currently not used
        void getLeaves_fast(
                size_t nodeId, size_t factId,
                PredId_t nodePred,
                std::shared_ptr<const TGSegment> data,
                size_t ruleIdx,
                size_t step,
                const std::vector<size_t> &incomingEdges,
                std::vector<std::vector<Literal>> &out);

        void exportEDBNode(JSON &out, Literal &l, size_t factId);

        void exportEDBNode(Literal &l, size_t factId, std::vector<Literal> &out);

        void getMappings(const Literal &l,
                const std::vector<uint64_t> &row,
                std::vector<std::pair<Term_t, Term_t>> &mappings);

        Literal ground(const Literal &l, const std::vector<
                std::pair<Term_t, Term_t>> &mappings,
                bool &ok);

        JSON getDerivationTree(size_t nodeId, size_t factId,
                DuplicateChecker *checker);

        std::vector<Term_t> getLeavesInDerivationTree(
                size_t nodeId,
                size_t factId,
                std::vector<std::vector<Literal>> &out,
                bool &validProof,
                DuplicateChecker *checker);

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
                const std::vector<size_t> &incomingEdges,
                DuplicateChecker *checker);

        void getLeavesInDerivationTree(
                size_t nodeId,
                size_t factId,
                std::vector<std::vector<Literal>> &out);

        std::vector<std::string> getListPredicates() const;

        std::string getTermText(Term_t t) const;

        JSON getNodeDetailsWithPredicate(std::string predName) const;

        JSON getNodeFacts(size_t nodeId) const;

        std::map<std::string, std::vector<std::vector<std::string>>> getAllFacts() const;

        std::pair<std::vector<std::pair<size_t, size_t>>, std::vector<Term_t>>
            getAllFactsPredicate(std::string predName) const;

        bool checkSoundnessDerivationTree(JSON &root, size_t threshold = ~0ul);

};

#endif
