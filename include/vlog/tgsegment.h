#ifndef _TG_SEGMENT_H
#define _TG_SEGMENT_H

#include <vlog/concepts.h>
#include <vlog/column.h>
#include <vlog/tgsegmentitr.h>
#include <vlog/tgsegmentcache.h>

#include <vector>
#include <map>


class TGSegment {
    private:

    public:
        virtual size_t getNRows() const = 0;

        virtual size_t getNColumns() const = 0;

        virtual bool isEmpty() const = 0;

        virtual std::unique_ptr<TGSegmentItr> iterator() const = 0;

        virtual std::unique_ptr<TGSegmentItr> sortBy(uint8_t field) const = 0;

        virtual std::unique_ptr<TGSegment> unique() const = 0;

        virtual std::unique_ptr<TGSegment> sort() const = 0;

        virtual std::unique_ptr<TGSegment> swap() const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendToWithProv(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out,
                size_t prov) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendToWithProv(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out,
                size_t prov) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual ~TGSegment() {}
};

class TGSegmentLegacy : public TGSegment {
    private:
        const size_t nrows;
        const bool sorted;
        const std::vector<std::shared_ptr<Column>> columns;

    public:
        TGSegmentLegacy(const std::vector<std::shared_ptr<Column>> &columns,
                size_t nrows, bool sorted=false) : nrows(nrows), sorted(sorted),
        columns(columns) {
        }

        size_t getNRows() const {
            return nrows;
        }

        size_t getNColumns() const {
            return columns.size();
        }

        bool isEmpty() const {
            return nrows == 0;
        }

        std::unique_ptr<TGSegmentItr> iterator() const;

        std::unique_ptr<TGSegmentItr> sortBy(uint8_t field) const;

        std::unique_ptr<TGSegment> unique() const;

        std::unique_ptr<TGSegment> sort() const;

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const;

        void appendToWithProv(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out,
                size_t prov) const;

        void appendTo(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const;

        void appendToWithProv(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<BinWithProv> &out,
                size_t prov) const;

        std::unique_ptr<TGSegment> swap() const;
};

template<typename K>
class TGSegmentImpl : public TGSegment {
    protected:
        std::vector<K> tuples;
        const size_t nodeId;
        const bool isSorted;

    public:
        TGSegmentImpl(std::vector<K> &t, const size_t nodeId, bool isSorted=false) :
            nodeId(nodeId), isSorted(isSorted) {
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

        virtual ~TGSegmentImpl() {}
};

template<typename S, typename K, typename I>
class UnaryTGSegmentImpl : public TGSegmentImpl<K> {
    public:
        UnaryTGSegmentImpl(std::vector<K> &tuples,
                const size_t nodeId) : TGSegmentImpl<K>(tuples, nodeId) {}

        size_t getNColumns() const {
            return 1;
        }

        std::unique_ptr<TGSegmentItr> iterator() const {
            return std::unique_ptr<TGSegmentItr>(new I(TGSegmentImpl<K>::getNodeId(), TGSegmentImpl<K>::tuples));
        }

        std::unique_ptr<TGSegmentItr> sortBy(uint8_t field) const {
            assert(field == 0);
            if (TGSegmentImpl<K>::isSorted) {
                return iterator();
            } else {
                SegmentCache &cache = SegmentCache::getInstance();
                const K* key = TGSegmentImpl<K>::tuples.data();
                if (!cache.contains(key)) {
                    std::vector<K> sortedTuples(TGSegmentImpl<K>::tuples);
                    std::sort(sortedTuples.begin(), sortedTuples.end());
                    cache.insert(key, sortedTuples);
                }
                const std::vector<K> &cachedSegment = cache.get(key);
                return std::unique_ptr<TGSegmentItr>(new I(TGSegmentImpl<K>::getNodeId(), cachedSegment));
            }
        }

        std::unique_ptr<TGSegment> unique() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            auto itr = std::unique(t.begin(), t.end());
            t.erase(itr, t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId()));
        }

        std::unique_ptr<TGSegment> sort() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            std::sort(t.begin(), t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId()));
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
                const size_t nodeId) : TGSegmentImpl<K>(tuples, nodeId) {}

        size_t getNColumns() const {
            return 2;
        }

        std::unique_ptr<TGSegmentItr> iterator() const {
            return std::unique_ptr<TGSegmentItr>(new I(TGSegmentImpl<K>::getNodeId(), TGSegmentImpl<K>::tuples));
        }

        std::unique_ptr<TGSegmentItr> sortBy(uint8_t field) const {
            if (TGSegmentImpl<K>::isSorted && field == 0) {
                return iterator();
            } else {
                SegmentCache &cache = SegmentCache::getInstance();
                const K* key = TGSegmentImpl<K>::tuples.data();
                if (!cache.contains(key, field)) {
                    std::vector<K> sortedTuples(TGSegmentImpl<K>::tuples);
                    if (field == 0) {
                        std::sort(sortedTuples.begin(), sortedTuples.end());
                    } else {
                        assert(field == 1);
                        std::sort(sortedTuples.begin(), sortedTuples.end(), invertedSorter<K>);
                    }
                    cache.insert(key, field, sortedTuples);
                }
                const std::vector<K> &cachedSegment = cache.get(key, field);
                return std::unique_ptr<TGSegmentItr>(new I(TGSegmentImpl<K>::getNodeId(), cachedSegment));
            }
        }

        std::unique_ptr<TGSegment> unique() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            auto itr = std::unique(t.begin(), t.end());
            t.erase(itr, t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId()));
        }

        std::unique_ptr<TGSegment> sort() const {
            auto t = std::vector<K>(TGSegmentImpl<K>::tuples);
            std::sort(t.begin(), t.end());
            return std::unique_ptr<TGSegment>(new S(t, TGSegmentImpl<K>::getNodeId()));
        }

};

class UnaryTGSegment : public UnaryTGSegmentImpl<UnaryTGSegment, Term_t, UnaryTGSegmentItr> {
    public:
        UnaryTGSegment(std::vector<Term_t> &tuples, const size_t nodeId) : UnaryTGSegmentImpl(tuples, nodeId) { }
        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            std::copy(tuples.begin(), tuples.end(), std::back_inserter(out));
        }

};

class UnaryWithProvTGSegment : public UnaryTGSegmentImpl<UnaryWithProvTGSegment, std::pair<Term_t,Term_t>, UnaryWithProvTGSegmentItr> {
    public:
        UnaryWithProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId) :
            UnaryTGSegmentImpl(tuples, nodeId) { }
        void appendToWithProv(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out,
                size_t prov) const {
            std::copy(tuples.begin(), tuples.end(), std::back_inserter(out));
        }
};

class UnaryWithConstProvTGSegment : public UnaryTGSegmentImpl<UnaryWithConstProvTGSegment, Term_t, UnaryTGSegmentItr> {
    public:
        UnaryWithConstProvTGSegment(std::vector<Term_t> &tuples,
                const size_t nodeId) :
            UnaryTGSegmentImpl(tuples, nodeId) { }
};

class BinaryTGSegment : public BinaryTGSegmentImpl<BinaryTGSegment, std::pair<Term_t,Term_t>,BinaryTGSegmentItr> {
    public:
        BinaryTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId) :
            BinaryTGSegmentImpl(tuples, nodeId) { }
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
    public:
        BinaryWithProvTGSegment(std::vector<BinWithProv> &tuples,
                const size_t nodeId) :
            BinaryTGSegmentImpl(tuples, nodeId) { }
        void appendToWithProv(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out,
                size_t prov) const {
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
};

class BinaryWithConstProvTGSegment : public BinaryTGSegmentImpl<BinaryWithConstProvTGSegment, std::pair<Term_t,Term_t>, BinaryTGSegmentItr> {
    public:
        BinaryWithConstProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId) :
            BinaryTGSegmentImpl(tuples, nodeId) { }
};

#endif
