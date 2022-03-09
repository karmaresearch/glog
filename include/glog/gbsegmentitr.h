#ifndef _TG_SEGMENT_ITR_H
#define _TG_SEGMENT_ITR_H

#include <vlog/column.h>

typedef enum SegProvenanceType {
    SEG_NOPROV, //no provenance
    SEG_SAMENODE, //all tuples come from the same node
    SEG_DIFFNODES, //tuples from different nodes
    SEG_FULLPROV
} SegProvenanceType;

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

    BinWithOff() : first(0), second(0), off(0) {
    }

    BinWithOff(size_t first, size_t second) : first(first), second(second) {
    }

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

        virtual size_t getProvenanceOffset(size_t proofNr, int pos) const = 0;

        //virtual size_t getOffset() const = 0;

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

        virtual size_t getNProofs() const
        {
            return 1;
        }

        virtual ~TGSegmentItr() {}
};

class TGSegmentCompProvItr : public TGSegmentItr {
    private:
        const SegProvenanceType provenanceType;
        const size_t nodeId;
        const size_t arity;
        const std::vector<Term_t> &data;

        int64_t m_currentIdx;
        int64_t currentIdx;

        bool shouldTrackProvenance() const {
            return provenanceType != SEG_NOPROV;
        }

    public:
        TGSegmentCompProvItr(const SegProvenanceType provenanceType,
                size_t nodeId, size_t arity,
                const std::vector<Term_t> &data) :
            provenanceType(provenanceType), nodeId(nodeId),
            arity(arity), data(data)
    {
        currentIdx = -arity;
    }

        bool hasNext() {
            return currentIdx < (int64_t)(data.size() - arity);
        }

        void next() {
            currentIdx += arity;
        }

        void mark() {
            m_currentIdx = currentIdx;
        }

        void reset() {
            currentIdx = m_currentIdx;
        }

        Term_t get(const int colIdx) {
            return data[currentIdx + colIdx];
        }

        size_t getNodeId() const {
            return nodeId;
        }

        size_t getProvenanceOffset(size_t proofNr, int pos) const {
            throw 10;
        }

        int getNFields() const {
            return arity;
        }

        size_t getNProofs() const {
            throw 10; //to implement
        }
};

class TGSegmentLegacyItr : public TGSegmentItr {
    private:
        const SegProvenanceType provenanceType;
        const std::vector<std::shared_ptr<Column>> &columns;
        const size_t nprovcolumns;
        std::vector<std::shared_ptr<ColumnReader>> readers;
        std::vector<Term_t> values;
        std::vector<Term_t> m_values;
        //int64_t counter;
        //int64_t m_counter;

        bool shouldTrackProvenance() const {
            return provenanceType != SEG_NOPROV;
        }

    public:
        TGSegmentLegacyItr(const std::vector<std::shared_ptr<Column>> &columns,
                const SegProvenanceType provenanceType,
                size_t nprovcolumns) :
            provenanceType(provenanceType),
            columns(columns),
            nprovcolumns(nprovcolumns)//, counter(-1)
    {
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
            //counter++;
        }

        void mark() {
            for(size_t i = 0; i < readers.size(); ++i) {
                readers[i]->mark();
                m_values[i] = values[i];
            }
            //m_counter = counter;
        }

        void reset() {
            for(size_t i = 0; i < readers.size(); ++i) {
                readers[i]->reset();
                values[i] = m_values[i];
            }
            //counter = m_counter;
        }

        Term_t get(const int colIdx) {
            return values[colIdx];
        }

        size_t getNodeId() const {
            return values[columns.size() - nprovcolumns];
        }

        /*size_t getOffset() const {
          return counter;
          }*/

        size_t getProvenanceOffset(size_t proofNr, int pos) const {
            assert(proofNr == 0);
            return values[columns.size() - nprovcolumns + pos + 1];
        }

        int getNFields() const {
            return shouldTrackProvenance() ? columns.size() - nprovcolumns : columns.size();
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

        //size_t getOffset() const {
        //    return currentIdx;
        //}

        virtual size_t getProvenanceOffset(size_t proofNr, int pos) const {
            throw 10;
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

        size_t getProvenanceOffset(size_t proofNr, int pos) const {
            assert(proofNr == 0);
            if (pos == 0)
                return tuples[currentIdx].second;
            else
                throw 10;
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

        size_t getProvenanceOffset(size_t proofNr, int pos) const {
            assert(proofNr == 0);
            if (pos == 0)
                return tuples[currentIdx].prov;
            else
                throw 10;
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

        size_t getProvenanceOffset(size_t proofNr, int pos) const {
            assert(proofNr == 0);
            if (pos == 0)
                return tuples[currentIdx].off;
            else
                throw 10;
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

        size_t getProvenanceOffset(size_t proofNr, int pos) const {
            assert(proofNr == 0);
            if (pos == 0)
                return tuples[currentIdx].prov;
            else
                throw 10;
        }

};

#endif
