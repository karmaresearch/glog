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

        virtual std::unique_ptr<TGSegmentDirectItr> directIterator() const = 0;

        virtual bool isSortedBy(std::vector<uint8_t> &fields) const = 0;

        virtual std::shared_ptr<TGSegment> sort() const = 0;

        virtual std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const = 0;

        virtual std::shared_ptr<TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const = 0;

        virtual std::shared_ptr<TGSegment> unique() const = 0;

        virtual std::shared_ptr<TGSegment> slice(const size_t nodeId,
                const size_t start,
                const size_t end) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual std::shared_ptr<TGSegment> shuffle(
                const std::vector<size_t> &idxs) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual std::string getName() const {
            return "TGSegment";
        }

        virtual std::shared_ptr<TGSegment> swap() const {
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

        std::shared_ptr<TGSegment> slice(const size_t nodeId,
                const size_t start,
                const size_t end) const;

        std::unique_ptr<TGSegmentItr> iterator() const;

        std::unique_ptr<TGSegmentDirectItr> directIterator() const;

        bool isSortedBy(std::vector<uint8_t> &fields) const;

        std::shared_ptr<TGSegment> sort() const;

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const;

        std::shared_ptr<TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const;

        std::shared_ptr<TGSegment> unique() const;

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const;

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const;

        void appendTo(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const;

        void appendTo(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<BinWithProv> &out) const;

        std::shared_ptr<TGSegment> swap() const;
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

template<typename S, typename K, typename I>
class TGSegmentImpl : public TGSegment {
    protected:
        std::shared_ptr<std::vector<K>> tuples;
        const size_t nodeId;
        const bool isSorted;
        const uint8_t sortedField;

    public:
        TGSegmentImpl(std::vector<K> &t, const size_t nodeId, bool isSorted=false,
                uint8_t sortedField=0) :
            nodeId(nodeId), isSorted(isSorted), sortedField(sortedField) {
                tuples = std::shared_ptr<std::vector<K>>(new std::vector<K>());
                tuples->swap(t);
            }
        TGSegmentImpl(std::shared_ptr<std::vector<K>> t, const size_t nodeId, bool isSorted=false,
                uint8_t sortedField=0) :
            tuples(t), nodeId(nodeId), isSorted(isSorted), sortedField(sortedField) {
            }

        size_t getNRows() const {
            return tuples->size();
        }

        const std::vector<K> &getTuples() const {
            return *tuples.get();
        }

        const size_t getNodeId() const {
            return nodeId;
        }

        bool isEmpty() const {
            return tuples->empty();
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

        virtual std::shared_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) const {
            if (start == 0 && end == getNRows()) {
                return std::shared_ptr<TGSegment>(new S(tuples,
                            nodeId,
                            TGSegmentImpl<S,K,I>::isSorted,
                            TGSegmentImpl<S,K,I>::sortedField));
            } else {
                std::vector<K> out(end - start);
                std::copy(TGSegmentImpl<S,K,I>::tuples->begin() + start, TGSegmentImpl<S,K,I>::tuples->begin() + end, out.begin());
                return std::shared_ptr<TGSegment>(new S(out,
                            nodeId,
                            TGSegmentImpl<S,K,I>::isSorted,
                            TGSegmentImpl<S,K,I>::sortedField));
            }
        }

        std::shared_ptr<TGSegment> shuffle(const std::vector<size_t> &idxs) const {
            std::vector<K> out;
            for(auto idx : idxs) out.push_back(TGSegmentImpl<S,K,I>::tuples->at(idx));
            return std::shared_ptr<TGSegment>(new S(out, TGSegmentImpl<S,K,I>::getNodeId(), false, 0));
        }

        std::unique_ptr<TGSegmentItr> iterator() const {
            return std::unique_ptr<TGSegmentItr>(new I(TGSegmentImpl<S,K,I>::getNodeId(), *TGSegmentImpl<S,K,I>::tuples.get()));
        }

        std::unique_ptr<TGSegmentDirectItr> directIterator() const {
            return std::unique_ptr<TGSegmentDirectItr>(new I(TGSegmentImpl<S,K,I>::getNodeId(), *TGSegmentImpl<S,K,I>::tuples.get()));
        }

        bool isSortedBy(std::vector<uint8_t> &fields) const {
            if (fields.size() != 1)
                return false;
            auto field = fields[0];
            assert(field == 0 || field == 1);
            return TGSegmentImpl<S,K,I>::isSorted && field == TGSegmentImpl<S,K,I>::sortedField;
        }

        std::shared_ptr<TGSegment> sort() const {
            auto t = std::vector<K>(*TGSegmentImpl<S,K,I>::tuples.get());
            std::sort(t.begin(), t.end());
            return std::shared_ptr<TGSegment>(new S(t, TGSegmentImpl<S,K,I>::getNodeId(), true, 0));
        }

        virtual std::shared_ptr<TGSegment> unique() const {
            auto t = std::vector<K>(*TGSegmentImpl<S,K,I>::tuples.get());
            auto itr = std::unique(t.begin(), t.end());
            t.erase(itr, t.end());
            return std::shared_ptr<TGSegment>(new S(t, TGSegmentImpl<S,K,I>::getNodeId(), true, 0));
        }

        virtual ~TGSegmentImpl() {}
};

template<typename S, typename K, typename I>
class UnaryTGSegmentImpl : public TGSegmentImpl<S,K,I> {
    public:
        UnaryTGSegmentImpl(std::vector<K> &tuples,
                const size_t nodeId, bool isSorted=false, uint8_t sortedField = 0)
            : TGSegmentImpl<S,K,I>(tuples, nodeId, isSorted, sortedField) {}

        UnaryTGSegmentImpl(std::shared_ptr<std::vector<K>> tuples,
                const size_t nodeId, bool isSorted=false, uint8_t sortedField = 0)
            : TGSegmentImpl<S,K,I>(tuples, nodeId, isSorted, sortedField) {}

        size_t getNColumns() const {
            return 1;
        }

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const {
            assert(fields.size() == 0 || (fields.size() == 1 && fields[0] == 0));
            std::vector<K> sortedTuples(*TGSegmentImpl<S,K,I>::tuples.get());
            std::sort(sortedTuples.begin(), sortedTuples.end());
            return std::shared_ptr<TGSegment>(
                    new S(sortedTuples, TGSegmentImpl<S,K,I>::getNodeId(), true, 0));
        }

        virtual std::string getName() const {
            return "TGUnarySegment";
        }
};

template<typename K>
bool invertedSorter(const K &a, const K &b) {
    return a.second < b.second || (a.second == b.second && a.first < b.first);
}

template<typename S, typename K, typename I>
class BinaryTGSegmentImpl : public TGSegmentImpl<S,K,I> {
    public:
        BinaryTGSegmentImpl(std::vector<K> &tuples,
                const size_t nodeId,
                bool isSorted=false, uint8_t sortedField = 0) : TGSegmentImpl<S,K,I>(tuples, nodeId,
                    isSorted, sortedField) {}

        BinaryTGSegmentImpl(std::shared_ptr<std::vector<K>> tuples,
                const size_t nodeId, bool isSorted=false, uint8_t sortedField = 0)
            : TGSegmentImpl<S,K,I>(tuples, nodeId, isSorted, sortedField) {}

        size_t getNColumns() const {
            return 2;
        }

        virtual std::string getName() const {
            return "TGBinarySegment";
        }

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const {
            uint8_t field = (fields.size() == 0 || fields[0] == 0) ? 0 : 1;
            std::vector<K> sortedTuples(*TGSegmentImpl<S,K,I>::tuples.get());
            if (field == 0) {
                std::sort(sortedTuples.begin(), sortedTuples.end());
            } else {
                assert(field == 1);
                std::sort(sortedTuples.begin(), sortedTuples.end(), invertedSorter<K>);
            }
            return std::shared_ptr<TGSegment>(
                    new S(sortedTuples, TGSegmentImpl<S,K,I>::getNodeId(), true, field));
        }

        std::shared_ptr<TGSegment> swap() const {
            std::vector<K> newtuples;
            for(const K &t : *TGSegmentImpl<S,K,I>::tuples.get()) {
                K newt = t;
                newt.first = t.second;
                newt.second = t.first;
                newtuples.push_back(newt);
            }
            return std::shared_ptr<TGSegment>(new S(newtuples, TGSegmentImpl<S,K,I>::getNodeId(), false, 0));
        }

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            if (colPos == 0) {
                for(const K &t : *TGSegmentImpl<S,K,I>::tuples.get()) {
                    out.push_back(t.first);
                }
            } else {
                assert(colPos == 1);
                for(const K &t : *TGSegmentImpl<S,K,I>::tuples.get()) {
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

        UnaryTGSegment(std::shared_ptr<std::vector<Term_t>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) {}

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
        }
};

class UnaryWithConstProvTGSegment : public UnaryTGSegmentImpl<UnaryWithConstProvTGSegment, Term_t, UnaryTGSegmentItr> {
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

            UnaryWithProvTGSegment(std::shared_ptr<std::vector<std::pair<Term_t,Term_t>>> tuples,
                    const size_t nodeId, bool isSorted, uint8_t sortedField) :
                UnaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

            void appendTo(uint8_t colPos,
                    std::vector<std::pair<Term_t, Term_t>> &out) const {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            }
            std::shared_ptr<TGSegment> unique() const {
                auto t = std::vector<std::pair<Term_t,Term_t>>(
                        *TGSegmentImpl<UnaryWithProvTGSegment, std::pair<Term_t,Term_t>, UnaryWithProvTGSegmentItr>::tuples.get());
                auto itr = std::unique(t.begin(), t.end(), cmpFirstTerm);
                t.erase(itr, t.end());
                return std::shared_ptr<TGSegment>(
                        new UnaryWithProvTGSegment(
                            t,
                            TGSegmentImpl<UnaryWithProvTGSegment,std::pair<Term_t,Term_t>,
                            UnaryWithProvTGSegmentItr>::getNodeId(),
                            true, 0));
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
                            isSorted,
                            sortedField));
            }
    };

class BinaryTGSegment : public BinaryTGSegmentImpl<BinaryTGSegment, std::pair<Term_t,Term_t>,BinaryTGSegmentItr> {
    public:
        BinaryTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        BinaryTGSegment(std::shared_ptr<std::vector<std::pair<Term_t,Term_t>>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else {
                assert(colPos1 == 1 && colPos2 == 0);
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.first));
                }
            }
        }
};

class BinaryWithConstProvTGSegment : public BinaryTGSegmentImpl<BinaryWithConstProvTGSegment, std::pair<Term_t,Term_t>, BinaryTGSegmentItr> {
    public:
        BinaryWithConstProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        BinaryWithConstProvTGSegment(std::shared_ptr<std::vector<std::pair<Term_t,Term_t>>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.first;
                    p.second = t.second;
                    p.node = nodeId;
                    out.push_back(p);
                }
            } else {
                assert(colPos1 == 1 && colPos2 == 0);
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.second;
                    p.second = t.first;
                    p.node = nodeId;
                    out.push_back(p);
                }
            }
        }

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            if (colPos == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.first, nodeId));
                }
            } else if (colPos == 1) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, nodeId));
                }
            } else {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            }
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            if (colPos1 != 0 || colPos2 != 1) {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            }
            std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
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

        BinaryWithProvTGSegment(std::shared_ptr<std::vector<BinWithProv>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else {
                assert(colPos1 == 1 && colPos2 == 0);
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.second;
                    p.second = t.first;
                    p.node = t.node;
                    out.push_back(p);
                }
            }
        }

        std::shared_ptr<TGSegment> unique() const {
            auto t = std::vector<BinWithProv>(
                    *TGSegmentImpl<BinaryWithProvTGSegment, BinWithProv, BinaryWithProvTGSegmentItr>::tuples.get());
            auto itr = std::unique(t.begin(), t.end(), cmpFirstSecondTerm);
            t.erase(itr, t.end());
            return std::shared_ptr<TGSegment>(new BinaryWithProvTGSegment(t,
                        TGSegmentImpl<
                        BinaryWithProvTGSegment,
                        BinWithProv,
                        BinaryWithProvTGSegmentItr>::getNodeId(), true, 0));
        }

        std::shared_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) const {
            std::vector<std::pair<Term_t,Term_t>> out(end - start);
            size_t m = 0;
            for(size_t j = start; j < end; ++j) {
                out[m].first = tuples->at(j).first;
                out[m].second = tuples->at(j).second;
                m++;
            }
            return std::shared_ptr<TGSegment>(
                    new BinaryWithConstProvTGSegment(out,
                        nodeId,
                        isSorted,
                        sortedField));
        }
};

#endif
