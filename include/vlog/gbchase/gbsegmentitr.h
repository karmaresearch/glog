#ifndef _TG_SEGMENT_ITR_H
#define _TG_SEGMENT_ITR_H

#include <vlog/column.h>

struct BinWithProv {
    size_t first, second, node;

    bool operator <(const BinWithProv &rhs) const {
        return first < rhs.first || (first == rhs.first && second < rhs.second) ||
            (first == rhs.first && second == rhs.second && node < rhs.node);
    }

    bool operator==(const BinWithProv& rhs) const {
        return first == rhs.first && second == rhs.second && node == rhs.node;
    }
};

struct BinWithOff {
    size_t first, second, off;

    bool operator <(const BinWithOff &rhs) const {
        return first < rhs.first || (first == rhs.first && second < rhs.second) ||
            (first == rhs.first && second == rhs.second && off < rhs.off);
    }

    bool operator==(const BinWithOff & rhs) const {
        return first == rhs.first && second == rhs.second && off == rhs.off;
    }
};


struct UnWithFullProv {
    size_t first, node, prov;

    bool operator <(const UnWithFullProv &rhs) const {
        return first < rhs.first ||
            (first == rhs.first && node < rhs.node) ||
            (first == rhs.first && node == rhs.node && prov < rhs.prov);
    }

    bool operator==(const UnWithFullProv& rhs) const {
        return first == rhs.first && node == rhs.node && prov == rhs.prov;
    }
};

struct BinWithFullProv {
    size_t first, second, node, prov;

    bool operator <(const BinWithFullProv &rhs) const {
        return first < rhs.first || (first == rhs.first && second < rhs.second) ||
            (first == rhs.first && second == rhs.second && node < rhs.node) ||
            (first == rhs.first && second == rhs.second && node == rhs.node && prov < rhs.prov);
    }

    bool operator==(const BinWithFullProv& rhs) const {
        return first == rhs.first && second == rhs.second && node == rhs.node
            && prov == rhs.prov;
    }
};

class TGSegmentItr {
    public:
        virtual bool hasNext() = 0;

        virtual void next() = 0;

        virtual Term_t get(const int colIdx) = 0;

        virtual size_t getNodeId() const = 0;

        virtual size_t getProvenanceOffset(int pos) const = 0;

        virtual int getNFields() const = 0;

        virtual void mark() = 0;

        virtual void reset() = 0;

        static int cmp(TGSegmentItr* inputLeft,
                TGSegmentItr* inputRight) {
            auto n = inputLeft->getNFields();
            for(size_t i = 0; i < n; ++i) {
                auto vl = inputLeft->get(i);
                auto vr = inputRight->get(i);
                if (vl < vr) {
                    return -1;
                } else if (vl > vr) {
                    return 1;
                }
            }
            return 0;
        }

        static int cmp(TGSegmentItr *inputLeft,
                TGSegmentItr *inputRight,
                std::vector<std::pair<int, int>> &joinVarPos) {
            for (int i = 0; i < joinVarPos.size(); i++) {
                const auto valLeft = inputLeft->get(joinVarPos[i].first);
                const auto valRight = inputRight->get(joinVarPos[i].second);
                if (valLeft < valRight)
                    return -1;
                else if (valLeft > valRight)
                    return 1;
            }
            return 0;
        }

        virtual ~TGSegmentItr() {}
};

class TGSegmentLegacyItr : public TGSegmentItr {
    private:
        const bool trackProvenance;
        const std::vector<std::shared_ptr<Column>> &columns;
        std::vector<std::shared_ptr<ColumnReader>> readers;
        std::vector<Term_t> values;
        std::vector<Term_t> m_values;

    public:
        TGSegmentLegacyItr(const std::vector<std::shared_ptr<Column>> &columns,
                const bool trackProvenance) :
            trackProvenance(trackProvenance),
            columns(columns) {
                for (const auto &c : columns)
                    readers.push_back(c->getReader());
                values.resize(columns.size());
                m_values.resize(columns.size());
            }

        bool hasNext() {
            bool o = true;
            for(const auto &r : readers) {
                o = r->hasNext();
                if (!o)
                    break;
            }
            return o;
        }

        void next() {
            for(size_t i = 0; i < readers.size(); ++i) {
                values[i] = readers[i]->next();
            }
        }

        void mark() {
            for(size_t i = 0; i < readers.size(); ++i) {
                readers[i]->mark();
                m_values[i] = values[i];
            }
        }

        void reset() {
            for(size_t i = 0; i < readers.size(); ++i) {
                readers[i]->reset();
                values[i] = m_values[i];
            }
        }

        Term_t get(const int colIdx) {
            return values[colIdx];
        }

        size_t getNodeId() const {
            return values.back();
        }

        size_t getProvenanceOffset(int pos) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }


        int getNFields() const {
            return trackProvenance ? columns.size() - 1 : columns.size();
        }
};

template<typename K>
class TGSegmentImplItr : public TGSegmentItr {
    private:
        const size_t nodeId;
        const int64_t nrows;
        int64_t m_currentIdx;
        std::shared_ptr<const TGSegment> selfref;

    protected:
        const std::vector<K> &tuples;
        int64_t currentIdx;

    public:
        TGSegmentImplItr(const size_t nodeId,
                const std::vector<K> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            nodeId(nodeId), nrows(tuples.size()),
            currentIdx(-1),
            selfref(selfref),
            tuples(tuples) {}

        void mark() {
            m_currentIdx = currentIdx;
        }

        void reset() {
            currentIdx = m_currentIdx;
        }

        bool hasNext() {
            return currentIdx < (nrows - 1);
        }

        void next() {
            currentIdx++;
        }

        virtual size_t getNodeId() const {
            return nodeId;
        }

        virtual size_t getProvenanceOffset(int pos) const {
            return currentIdx;
        }

        virtual ~TGSegmentImplItr() {}
};

class UnaryTGSegmentItr : public TGSegmentImplItr<Term_t> {
    public:
        UnaryTGSegmentItr(const size_t nodeId,
                const std::vector<Term_t> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            return tuples[currentIdx];
        }

        int getNFields() const { return 1; }
};

class UnaryWithProvTGSegmentItr : public TGSegmentImplItr<std::pair<Term_t,Term_t>> {
    public:
        UnaryWithProvTGSegmentItr(const size_t nodeId,
                const std::vector<std::pair<Term_t,Term_t>> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            return tuples[currentIdx].first;
        }

        int getNFields() const { return 1; }

        size_t getNodeId() const {
            return tuples[currentIdx].second;
        }
};

class UnaryWithConstNodeFullProvTGSegmentItr : public TGSegmentImplItr<std::pair<Term_t,Term_t>> {
    public:
        UnaryWithConstNodeFullProvTGSegmentItr(const size_t nodeId,
                const std::vector<std::pair<Term_t,Term_t>> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            return tuples[currentIdx].first;
        }

        int getNFields() const { return 1; }

        size_t getProvenanceOffset() const {
            return tuples[currentIdx].second;
        }
};

class UnaryWithFullProvTGSegmentItr : public TGSegmentImplItr<UnWithFullProv> {
    public:
        UnaryWithFullProvTGSegmentItr(const size_t nodeId,
                const std::vector<UnWithFullProv> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            return tuples[currentIdx].first;
        }

        int getNFields() const { return 1; }

        size_t getNodeId() const {
            return tuples[currentIdx].node;
        }

        size_t getProvenanceOffset(int pos) const {
            assert(pos == 0);
            return tuples[currentIdx].prov;
        }
};

class BinaryTGSegmentItr : public TGSegmentImplItr<std::pair<Term_t,Term_t>> {
    public:
        BinaryTGSegmentItr(const size_t nodeId,
                const std::vector<std::pair<Term_t,Term_t>> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            if (colIdx == 0) {
                return tuples[currentIdx].first;
            } else {
                return tuples[currentIdx].second;
            }
        }

        int getNFields() const { return 2; }
};

class BinaryWithProvTGSegmentItr : public TGSegmentImplItr<BinWithProv> {
    public:
        BinaryWithProvTGSegmentItr(const size_t nodeId,
                const std::vector<BinWithProv> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            if (colIdx == 0) {
                return tuples[currentIdx].first;
            } else {
                return tuples[currentIdx].second;
            }
        }

        int getNFields() const { return 2; }

        size_t getNodeId() const {
            return tuples[currentIdx].node;
        }
};

class BinaryWithOffTGSegmentItr : public TGSegmentImplItr<BinWithOff> {
    public:
        BinaryWithOffTGSegmentItr(const size_t nodeId,
                const std::vector<BinWithOff> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            if (colIdx == 0) {
                return tuples[currentIdx].first;
            } else {
                return tuples[currentIdx].second;
            }
        }

        int getNFields() const { return 2; }

        size_t getProvenanceOffset() const {
            return tuples[currentIdx].off;
        }
};


class BinaryWithFullProvTGSegmentItr : public TGSegmentImplItr<BinWithFullProv> {
    public:
        BinaryWithFullProvTGSegmentItr(const size_t nodeId,
                const std::vector<BinWithFullProv> &tuples,
                std::shared_ptr<const TGSegment> selfref = NULL) :
            TGSegmentImplItr(nodeId, tuples, selfref) {}

        Term_t get(const int colIdx) {
            if (colIdx == 0) {
                return tuples[currentIdx].first;
            } else {
                return tuples[currentIdx].second;
            }
        }

        int getNFields() const { return 2; }

        size_t getNodeId() const {
            return tuples[currentIdx].node;
        }

        size_t getProvenanceOffset(int pos) const {
            assert(pos == 0);
            return tuples[currentIdx].prov;
        }
};

#endif
