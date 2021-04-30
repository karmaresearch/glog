#ifndef _GBSEGMENT_UNARY_H
#define _GBSEGMENT_UNARY_H

#include <vlog/gbchase/gbsegmentitr.h>
#include <vlog/gbchase/gbsegment.h>

template<typename S, typename K, typename I, SegProvenanceType CP>
class UnaryTGSegmentImpl : public TGSegmentImpl<S,K,I,CP> {
    public:
        UnaryTGSegmentImpl(std::vector<K> &tuples,
                const size_t nodeId, bool isSorted=false, uint8_t sortedField = 0)
            : TGSegmentImpl<S,K,I,CP>(tuples, nodeId, isSorted, sortedField) {}

        UnaryTGSegmentImpl(std::shared_ptr<std::vector<K>> tuples,
                const size_t nodeId, bool isSorted=false, uint8_t sortedField = 0)
            : TGSegmentImpl<S,K,I,CP>(tuples, nodeId, isSorted, sortedField) {}

        size_t getNColumns() const {
            return 1;
        }

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const {
            assert(fields.size() == 0 || (fields.size() == 1 && fields[0] == 0));
            std::vector<K> sortedTuples(*TGSegmentImpl<S,K,I,CP>::tuples.get());
            std::sort(sortedTuples.begin(), sortedTuples.end());
            return std::shared_ptr<TGSegment>(
                    new S(sortedTuples, TGSegmentImpl<S,K,I,CP>::getNodeId(), true, 0));
        }

        virtual std::string getName() const {
            return "TGUnarySegment";
        }
};

class UnaryTGSegment : public UnaryTGSegmentImpl<UnaryTGSegment, Term_t, UnaryTGSegmentItr, SEG_NOPROV> {
    public:
        UnaryTGSegment(std::vector<Term_t> &tuples, const size_t nodeId,
                bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        UnaryTGSegment(std::shared_ptr<std::vector<Term_t>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) {}

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            assert(colPos1 == colPos2 == 0);
            for(const auto &t : *tuples.get()) {
                out.push_back(std::make_pair(t, t));
            }
        }
};

class UnaryWithConstProvTGSegment : public UnaryTGSegmentImpl<
                                    UnaryWithConstProvTGSegment,
                                    Term_t,
                                    UnaryTGSegmentItr,
                                    SEG_SAMENODE>
{
    public:
        UnaryWithConstProvTGSegment(std::vector<Term_t> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        UnaryWithConstProvTGSegment(std::shared_ptr<std::vector<Term_t>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) {}

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            for(const auto &value : *tuples.get())
                out.push_back(std::make_pair(value, nodeId));
        }

        void appendTo(uint8_t colPos,
                std::vector<Term_t> &out) const {
            std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            assert(colPos1 == colPos2);
            assert(colPos1 == 0);
            for(const auto &value : *tuples.get()) {
                BinWithProv p;
                p.first = value;
                p.second = value;
                p.node = nodeId;
                out.push_back(p);
            }
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            assert(colPos1 == colPos2);
            assert(colPos1 == 0);
            for(const auto &value : *tuples.get()) {
                out.push_back(std::make_pair(value, value));
            }
        }

        size_t countHits(const std::vector<Term_t> &terms,
                int column) const {
            assert(column == 0);
            size_t c = 0;
            for(auto &t : terms) {
                if (std::binary_search(tuples->begin(),
                            tuples->end(), t))
                    c++;
            }
            return c;
        }
};

class UnaryWithProvTGSegment : public UnaryTGSegmentImpl<UnaryWithProvTGSegment,
    std::pair<Term_t,Term_t>, UnaryWithProvTGSegmentItr, SEG_DIFFNODES>
{
    private:
        static bool cmpFirstTerm(const std::pair<Term_t,Term_t> &a,
                const std::pair<Term_t, Term_t> &b) {
            return a.first == b.first;
        }
        static bool sortSecondTerm(const std::pair<Term_t,Term_t> &a,
                const std::pair<Term_t, Term_t> &b) {
            return a.second < b.second || (a.second == b.second &&
                    a.first < b.first);
        }
    public:
        UnaryWithProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        UnaryWithProvTGSegment(std::shared_ptr<std::vector<std::pair<Term_t,Term_t>>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            assert(colPos1 == colPos2 == 0);
            for(const auto &value : *tuples.get()) {
                BinWithProv p;
                p.first = value.first;
                p.second = value.first;
                p.node = value.second;
                out.push_back(p);
            }
        }

        std::shared_ptr<const TGSegment> unique() const {
            auto t = std::vector<std::pair<Term_t,Term_t>>(
                    *TGSegmentImpl<UnaryWithProvTGSegment,
                    std::pair<Term_t,Term_t>, UnaryWithProvTGSegmentItr, SEG_DIFFNODES>::tuples.get());
            auto itr = std::unique(t.begin(), t.end(), cmpFirstTerm);
            t.erase(itr, t.end());
            return std::shared_ptr<const TGSegment>(
                    new UnaryWithProvTGSegment(
                        t,
                        TGSegmentImpl<UnaryWithProvTGSegment,
                        std::pair<Term_t, Term_t>,
                        UnaryWithProvTGSegmentItr, SEG_DIFFNODES>::getNodeId(),
                        true, 0));
        }

        std::shared_ptr<const TGSegment> sortByProv() const {
            auto t = std::vector<std::pair<Term_t,Term_t>>(
                    *TGSegmentImpl<UnaryWithProvTGSegment,
                    std::pair<Term_t,Term_t>,
                    UnaryWithProvTGSegmentItr, SEG_DIFFNODES>::tuples.get());
            std::sort(t.begin(), t.end(), sortSecondTerm);
            return std::shared_ptr<const TGSegment>(
                    new UnaryWithProvTGSegment(
                        t,
                        TGSegmentImpl<UnaryWithProvTGSegment,
                        std::pair<Term_t, Term_t>,
                        UnaryWithProvTGSegmentItr, SEG_DIFFNODES>::getNodeId(),
                        false, 0));
        }

        std::shared_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) const {
            std::vector<Term_t> out(end - start);
            size_t m = 0;
            for(size_t j = start; j < end; ++j) {
                out[m++] = tuples->at(j).first;
            }
            return std::shared_ptr<TGSegment>(
                    new UnaryWithConstProvTGSegment(out,
                        nodeId,
                        f_isSorted,
                        sortedField));
        }
};

/* Both the node and the offset can be determined automatically. The container stores only the constants */
class UnaryWithConstNodeOffFullProvTGSegment : public UnaryTGSegmentImpl<UnaryWithConstNodeOffFullProvTGSegment,
    Term_t, UnaryTGSegmentItr, SEG_FULLPROV>
{
    public:
        UnaryWithConstNodeOffFullProvTGSegment(std::vector<Term_t> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        UnaryWithConstNodeOffFullProvTGSegment(std::shared_ptr<std::vector<Term_t>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }


        std::shared_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) const {
            throw 10;
        }
};

/* The node is constant. The container stores the constants and the offset. */
class UnaryWithConstNodeFullProvTGSegment : public UnaryTGSegmentImpl<UnaryWithConstNodeFullProvTGSegment,
    std::pair<Term_t,Term_t>, UnaryWithConstNodeFullProvTGSegmentItr, SEG_FULLPROV>
{
    private:
        static bool cmpFirstTerm(const std::pair<Term_t,Term_t> &a,
                const std::pair<Term_t, Term_t> &b) {
            return a.first == b.first;
        }
        static bool sortSecondTerm(const std::pair<Term_t,Term_t> &a,
                const std::pair<Term_t, Term_t> &b) {
            return a.second < b.second || (a.second == b.second &&
                    a.first < b.first);
        }
    public:
        UnaryWithConstNodeFullProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        UnaryWithConstNodeFullProvTGSegment(std::shared_ptr<std::vector<std::pair<Term_t,Term_t>>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithFullProv> &out) const {
            assert(colPos1 == colPos2 == 0);
            for(const auto &value : *tuples.get()) {
                BinWithFullProv p;
                p.first = value.first;
                p.second = value.first;
                p.node = getNodeId();
                p.prov = value.second;
                out.push_back(p);
            }
        }

        std::shared_ptr<const TGSegment> unique() const {
            auto t = std::vector<std::pair<Term_t,Term_t>>(
                    *TGSegmentImpl<UnaryWithConstNodeFullProvTGSegment,
                    std::pair<Term_t,Term_t>, UnaryWithConstNodeFullProvTGSegmentItr, SEG_FULLPROV>::tuples.get());
            auto itr = std::unique(t.begin(), t.end(), cmpFirstTerm);
            t.erase(itr, t.end());
            return std::shared_ptr<const TGSegment>(
                    new UnaryWithConstNodeFullProvTGSegment(
                        t,
                        TGSegmentImpl<UnaryWithConstNodeFullProvTGSegment,
                        std::pair<Term_t, Term_t>,
                        UnaryWithConstNodeFullProvTGSegmentItr, SEG_FULLPROV>::getNodeId(),
                        true, 0));
        }
};

/* The container stores the values, the nodes, and the offsets */
class UnaryWithFullProvTGSegment : public UnaryTGSegmentImpl<UnaryWithFullProvTGSegment,
    UnWithFullProv, UnaryWithFullProvTGSegmentItr, SEG_FULLPROV>
{
    private:
        static bool cmpFirstTerm(const UnWithFullProv &a,
                const UnWithFullProv &b) {
            return a.first == b.first;
        }
        static bool sortByNode(const UnWithFullProv &a,
                const UnWithFullProv &b) {
            return a.node < b.node || (a.node == b.node &&
                    a.first < b.first);
        }

    public:
        UnaryWithFullProvTGSegment(std::vector<UnWithFullProv> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        UnaryWithFullProvTGSegment(std::shared_ptr<std::vector<UnWithFullProv>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        void appendTo(uint8_t colPos,
                std::vector<UnWithFullProv> &out) const {
            std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithFullProv> &out) const {
            assert(colPos1 == colPos2 == 0);
            for(const auto &value : *tuples.get()) {
                BinWithFullProv p;
                p.first = value.first;
                p.second = value.first;
                p.node = value.node;
                p.prov = value.prov;
                out.push_back(p);
            }
        }

        std::shared_ptr<const TGSegment> unique() const {
            auto t = std::vector<UnWithFullProv>(
                    *TGSegmentImpl<UnaryWithFullProvTGSegment,
                    UnWithFullProv, UnaryWithFullProvTGSegmentItr, SEG_FULLPROV>::tuples.get());
            auto itr = std::unique(t.begin(), t.end(), cmpFirstTerm);
            t.erase(itr, t.end());
            return std::shared_ptr<const TGSegment>(
                    new UnaryWithFullProvTGSegment(
                        t,
                        TGSegmentImpl<UnaryWithFullProvTGSegment,
                        UnWithFullProv,
                        UnaryWithFullProvTGSegmentItr, SEG_FULLPROV>::getNodeId(),
                        true, 0));
        }

        std::shared_ptr<const TGSegment> sortByProv() const {
            auto t = std::vector<UnWithFullProv>(
                    *TGSegmentImpl<UnaryWithFullProvTGSegment,
                    UnWithFullProv,
                    UnaryWithFullProvTGSegmentItr, SEG_FULLPROV>::tuples.get());
            std::sort(t.begin(), t.end(), sortByNode);
            return std::shared_ptr<const TGSegment>(
                    new UnaryWithFullProvTGSegment(
                        t,
                        TGSegmentImpl<UnaryWithFullProvTGSegment,
                        UnWithFullProv,
                        UnaryWithFullProvTGSegmentItr, SEG_FULLPROV>::getNodeId(),
                        false, 0));
        }

        std::shared_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) const {
            throw 10;
        }
};

#endif
