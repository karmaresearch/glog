#ifndef _TG_SEGMENT_H
#define _TG_SEGMENT_H

#include <vlog/concepts.h>
#include <vlog/column.h>
#include <vlog/tgsegmentitr.h>

#include <vector>
#include <map>


class TGSegment {
    private:

    public:
        virtual size_t getNRows() const = 0;

        virtual size_t getNColumns() const = 0;

        virtual bool isEmpty() const = 0;

        virtual std::unique_ptr<TGSegmentItr> iterator() const = 0;

        virtual bool isSortedBy(uint8_t field) const = 0;

        virtual std::unique_ptr<TGSegment> sort() const = 0;

        virtual std::shared_ptr<TGSegment> sortBy(uint8_t field) const = 0;

        virtual std::shared_ptr<TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const = 0;

        virtual std::unique_ptr<TGSegment> unique() const = 0;

        virtual std::unique_ptr<TGSegment> slice(const size_t nodeId,
                const size_t start,
                const size_t end) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual std::unique_ptr<TGSegment> shuffle(
                const std::vector<size_t> &idxs) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual std::string getName() const {
            return "TGSegment";
        }

        virtual std::unique_ptr<TGSegment> swap() const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual ~TGSegment() {}
};

class TGSegmentLegacy : public TGSegment {
    private:
        const size_t nrows;
        const bool isSorted;
        const uint8_t sortedField;
        const std::vector<std::shared_ptr<Column>> columns;
        const bool trackProvenance;

    public:
        TGSegmentLegacy(const std::vector<std::shared_ptr<Column>> &columns,
                size_t nrows, bool isSorted=false, uint8_t sortedField = 0,
                bool trackProvenance = false) :
            nrows(nrows),
            isSorted(isSorted),
            sortedField(sortedField),
            columns(columns),
            trackProvenance(trackProvenance) {
            }

        size_t getNRows() const {
            return nrows;
        }

        size_t getNColumns() const {
            return trackProvenance ? columns.size() - 1 : columns.size();
        }

        bool isEmpty() const {
            return nrows == 0;
        }

        virtual std::string getName() const {
            return "TGSegmentLegacy";
        }

        std::unique_ptr<TGSegmentItr> iterator() const;

        bool isSortedBy(uint8_t field) const;

        std::unique_ptr<TGSegment> sort() const;

        std::shared_ptr<TGSegment> sortBy(uint8_t field) const;

        std::shared_ptr<TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const;

        std::unique_ptr<TGSegment> unique() const;

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const;

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const;

        void appendTo(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const;

        void appendTo(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<BinWithProv> &out) const;

        std::unique_ptr<TGSegment> swap() const;
};

struct ProvSorter {
    const size_t *tuples;
    const int ncols;

    ProvSorter(const size_t *tuples, const int ncols) :
        tuples(tuples), ncols(ncols) {
        }

    bool operator ()(const size_t &a, const size_t &b) const {
        const size_t ba = a * ncols;
        const size_t bb = b * ncols;
        for(int i = 0; i < ncols; ++i) {
            if (tuples[ba] < tuples[bb])
                return true;
        }
        return false;
    }
};

template<typename K>
class TGSegmentImpl : public TGSegment {
    protected:
        const std::shared_ptr<std::vector<K>> tuples;
        const K* rawpointer;
        const size_t nodeId;
        const bool isSorted;
        const uint8_t sortedField;

    public:
        TGSegmentImpl(std::vector<K> &t, const size_t nodeId, bool isSorted=false,
                uint8_t sortedField=0) :
            nodeId(nodeId), isSorted(isSorted), sortedField(sortedField) {
                tuples.swap(t);
            }

        size_t getNRows() const {
            return tuples.size();
        }

        const std::vector<K> &getTuples() const {
            return tuples;
        }

        const size_t getNodeId() const {
            return nodeId;
        }

        bool isEmpty() const {
            return tuples.empty();
        }

        virtual std::shared_ptr<TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const {
            const auto nrows = nodes.size() / ncols;
            idxs.resize(nrows);
            for(size_t i = 0; i < nrows; ++i) idxs[i] = i;
            std::sort(idxs.begin(), idxs.end(), ProvSorter(nodes.data(), ncols));
            return shuffle(idxs);
        }

        virtual ~TGSegmentImpl() {}
};

template<typename S, typename K, typename I>
class UnaryTGSegmentImpl : public TGSegmentImpl<K> {
    public:
        UnaryTGSegmentImpl(std::vector<K> &tuples,
                const size_t nodeId, bool isSorted=false, uint8_t sortedField = 0)
            : TGSegmentImpl<K>(tuples, nodeId, isSorted, sortedField) {}

        size_t getNColumns() const {
            return 1;
        }

        std::unique_ptr<TGSegmentItr> iterator() const {
            return std::unique_ptr<TGSegmentItr>(new I(TGSegmentImpl<K>::getNodeId(), TGSegmentImpl<K>::tuples));
        }

        bool isSortedBy(uint8_t field) const {
            assert(field == 0);
            return TGSegmentImpl<K>::isSorted;
        }

        std::shared_ptr<TGSegment> sortBy(uint8_t field) const {
            assert(field == 0);
            std::vector<K> sortedTuples(TGSegmentImpl<K>::tuples);
            std::sort(sortedTuples.begin(), sortedTuples.end());
            return std::shared_ptr<TGSegment>(new S(sortedTuples, TGSegmentImpl<K>::getNodeId(), true, field));
        }

        virtual std::string getName() const {
            return "TGUnarySegment";
        }

        virtual std::unique_ptr<TGSegment> unique() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            auto itr = std::unique(t.begin(), t.end());
            t.erase(itr, t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId(), true, 0));
        }

        std::unique_ptr<TGSegment> sort() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            std::sort(t.begin(), t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId(), true, 0));
        }

        std::unique_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) {
            std::vector<K> out(end - start);
            std::copy(TGSegmentImpl<K>::tuples.begin() + start, TGSegmentImpl<K>::tuples.begin() + end, out.begin());
            return std::unique_ptr<TGSegment>(new S(out,
                        nodeId,
                        TGSegmentImpl<K>::isSorted,
                        TGSegmentImpl<K>::sortedField));
        }

        std::unique_ptr<TGSegment> shuffle(const std::vector<size_t> &idxs) {
            std::vector<K> out;
            for(auto idx : idxs) out.push_back(TGSegmentImpl<K>::tuples[idx]);
            return std::unique_ptr<TGSegment>(new S(out, TGSegmentImpl<K>::getNodeId(), false, 0));
        }
};

template<typename K>
bool invertedSorter(const K &a, const K &b) {
    return a.second < b.second || (a.second == b.second && a.first < b.first);
}

template<typename S, typename K, typename I>
class BinaryTGSegmentImpl : public TGSegmentImpl<K> {
    public:
        BinaryTGSegmentImpl(std::vector<K> &tuples,
                const size_t nodeId,
                bool isSorted=false, uint8_t sortedField = 0) : TGSegmentImpl<K>(tuples, nodeId,
                    isSorted, sortedField) {}

        size_t getNColumns() const {
            return 2;
        }

        std::unique_ptr<TGSegmentItr> iterator() const {
            return std::unique_ptr<TGSegmentItr>(new I(TGSegmentImpl<K>::getNodeId(), TGSegmentImpl<K>::tuples));
        }

        bool isSortedBy(uint8_t field) const {
            assert(field == 0 || field == 1);
            return TGSegmentImpl<K>::isSorted && field == TGSegmentImpl<K>::sortedField;
        }

        virtual std::string getName() const {
            return "TGBinarySegment";
        }

        std::shared_ptr<TGSegment> sortBy(uint8_t field) const {
            std::vector<K> sortedTuples(TGSegmentImpl<K>::tuples);
            if (field == 0) {
                std::sort(sortedTuples.begin(), sortedTuples.end());
            } else {
                assert(field == 1);
                std::sort(sortedTuples.begin(), sortedTuples.end(), invertedSorter<K>);
            }
            return std::shared_ptr<TGSegment>(new S(sortedTuples, TGSegmentImpl<K>::getNodeId(), true, field));

        }

        virtual std::unique_ptr<TGSegment> unique() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            auto itr = std::unique(t.begin(), t.end());
            t.erase(itr, t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId(), true, 0));
        }

        std::unique_ptr<TGSegment> sort() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            std::sort(t.begin(), t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId(), true, 0));
        }

        std::unique_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) {
            std::vector<K> out(end - start);
            std::copy(TGSegmentImpl<K>::tuples.begin() + start, TGSegmentImpl<K>::tuples.begin() + end, out.begin());
            return std::unique_ptr<TGSegment>(new S(out,
                        nodeId,
                        TGSegmentImpl<K>::isSorted,
                        TGSegmentImpl<K>::sortedField));
        }

        std::unique_ptr<TGSegment> shuffle(const std::vector<size_t> &idxs) {
            std::vector<K> out;
            for(auto idx : idxs) out.push_back(TGSegmentImpl<K>::tuples[idx]);
            return std::unique_ptr<TGSegment>(new S(out, TGSegmentImpl<K>::getNodeId(), false, 0));
        }

        std::unique_ptr<TGSegment> swap() const {
            std::vector<K> newtuples;
            for(const K &t : TGSegmentImpl<K>::tuples) {
                K newt = t;
                newt.first = t.second;
                newt.second = t.first;
                newtuples.push_back(newt);
            }
            return std::unique_ptr<TGSegment>(new S(newtuples, TGSegmentImpl<K>::getNodeId(), false, 0));
        }

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            if (colPos == 0) {
                for(const K &t : TGSegmentImpl<K>::tuples) {
                    out.push_back(t.first);
                }
            } else {
                assert(colPos == 1);
                for(const K &t : TGSegmentImpl<K>::tuples) {
                    out.push_back(t.second);
                }
            }
        }
};

class UnaryTGSegment : public UnaryTGSegmentImpl<UnaryTGSegment, Term_t, UnaryTGSegmentItr> {
    public:
        UnaryTGSegment(std::vector<Term_t> &tuples, const size_t nodeId,
                bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }
        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            std::copy(tuples.begin(), tuples.end(), std::back_inserter(out));
        }
};

class UnaryWithConstProvTGSegment : public UnaryTGSegmentImpl<UnaryWithConstProvTGSegment, Term_t, UnaryTGSegmentItr> {
    public:
        UnaryWithConstProvTGSegment(std::vector<Term_t> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }
        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            for(const auto &value : tuples)
                out.push_back(std::make_pair(value, nodeId));
        }
};

class UnaryWithProvTGSegment : public UnaryTGSegmentImpl<UnaryWithProvTGSegment,
    std::pair<Term_t,Term_t>, UnaryWithProvTGSegmentItr> {

        private:
            static bool cmpFirstTerm(const std::pair<Term_t,Term_t> &a,
                    const std::pair<Term_t, Term_t> &b) {
                return a.first == b.first;
            }
        public:
            UnaryWithProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                    const size_t nodeId, bool isSorted, uint8_t sortedField) :
                UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }
            void appendTo(uint8_t colPos,
                    std::vector<std::pair<Term_t, Term_t>> &out) const {
                std::copy(tuples.begin(), tuples.end(), std::back_inserter(out));
            }
            std::unique_ptr<TGSegment> unique() const {
                auto t = std::vector<std::pair<Term_t,Term_t>>(
                        TGSegmentImpl<std::pair<Term_t,Term_t>>::tuples);
                auto itr = std::unique(t.begin(), t.end(), cmpFirstTerm);
                t.erase(itr, t.end());
                return std::unique_ptr<TGSegment>(
                        new UnaryWithProvTGSegment(
                            t,
                            TGSegmentImpl<std::pair<Term_t,Term_t>>::getNodeId(),
                            true, 0));
            }
            std::unique_ptr<TGSegment> slice(size_t nodeId,
                    const size_t start, const size_t end) {
                std::vector<Term_t> out(end - start);
                size_t m = 0;
                for(size_t j = start; j < end; ++j) {
                    out[m++] = tuples[j].first;
                }
                return std::unique_ptr<TGSegment>(
                        new UnaryWithConstProvTGSegment(out,
                            nodeId,
                            isSorted,
                            sortedField));
            }
    };


class BinaryTGSegment : public BinaryTGSegmentImpl<BinaryTGSegment, std::pair<Term_t,Term_t>,BinaryTGSegmentItr> {
    public:
        BinaryTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }
        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples.begin(), tuples.end(), std::back_inserter(out));
            } else {
                assert(colPos1 == 1 && colPos2 == 0);
                for(auto &t : tuples) {
                    out.push_back(std::make_pair(t.second, t.first));
                }
            }
        }
};

class BinaryWithProvTGSegment : public BinaryTGSegmentImpl<BinaryWithProvTGSegment, BinWithProv, BinaryWithProvTGSegmentItr> {
    private:
        static bool cmpFirstSecondTerm(const BinWithProv &a,
                const BinWithProv &b) {
            return a.first == b.first && a.second == b.second;
        }

    public:
        BinaryWithProvTGSegment(std::vector<BinWithProv> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }
        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples.begin(), tuples.end(), std::back_inserter(out));
            } else {
                assert(colPos1 == 1 && colPos2 == 0);
                for(auto &t : tuples) {
                    BinWithProv p;
                    p.first = t.second;
                    p.second = t.first;
                    p.node = t.node;
                    out.push_back(p);
                }
            }
        }
        std::unique_ptr<TGSegment> unique() const {
            auto t = std::vector<BinWithProv>(TGSegmentImpl<BinWithProv>::tuples);
            auto itr = std::unique(t.begin(), t.end(), cmpFirstSecondTerm);
            t.erase(itr, t.end());
            return std::unique_ptr<TGSegment>(new BinaryWithProvTGSegment(t, TGSegmentImpl<BinWithProv>::getNodeId(), true, 0));
        }
};

class BinaryWithConstProvTGSegment : public BinaryTGSegmentImpl<BinaryWithConstProvTGSegment, std::pair<Term_t,Term_t>, BinaryTGSegmentItr> {
    public:
        BinaryWithConstProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }
};

#endif
