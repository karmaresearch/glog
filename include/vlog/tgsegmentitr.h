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

class TGSegmentItr {
    public:
        virtual void mark() = 0;

        virtual void reset() = 0;

        virtual bool hasNext() = 0;

        virtual void next() = 0;

        virtual Term_t get(const int colIdx) = 0;

        virtual size_t getNodeId() const = 0;

        virtual int getNFields() const = 0;

        virtual ~TGSegmentItr() {}
};

class TGSegmentLegacyItr : public TGSegmentItr {
    private:
        const bool trackProvenance;
        const std::vector<std::shared_ptr<Column>> &columns;

        int64_t nrows;
        std::vector<const Term_t*> vectors;
        int64_t currentRowIdx;

        int64_t m_currentRowIdx;


    public:
        TGSegmentLegacyItr(const std::vector<std::shared_ptr<Column>> &columns,
                const bool trackProvenance) :
                trackProvenance(trackProvenance),
                columns(columns) {
            currentRowIdx = -1;
            if (columns.size() > 0) {
                nrows = columns[0]->size();
            } else {
                nrows = 1;
            }
            auto ncols = trackProvenance ? columns.size() - 1 : columns.size();
            for(int i = 0; i < ncols; ++i) {
                vectors.push_back(columns[i]->getVectorRef().data());
            }
        }

        void mark() {
            m_currentRowIdx = currentRowIdx;
        }

        void reset() {
            currentRowIdx = m_currentRowIdx;
        }

        bool hasNext() {
            return currentRowIdx < (nrows - 1);
        }

        void next() {
            currentRowIdx++;
        }

        Term_t get(const int colIdx) {
            return vectors[colIdx][currentRowIdx];
        }

        size_t getNodeId() const {
            return columns.back()->getValue(currentRowIdx);
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

    protected:
        const std::vector<K> &tuples;
        int64_t currentIdx;

    public:
        TGSegmentImplItr(const size_t nodeId,
                const std::vector<K> &tuples) : nodeId(nodeId), nrows(tuples.size()),
        currentIdx(-1),
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

        virtual ~TGSegmentImplItr() {}
};

class UnaryTGSegmentItr : public TGSegmentImplItr<Term_t> {
    public:
        UnaryTGSegmentItr(const size_t nodeId,
                const std::vector<Term_t> &tuples) :
            TGSegmentImplItr(nodeId, tuples) {}

        Term_t get(const int colIdx) {
            return tuples[currentIdx];
        }

        int getNFields() const { return 1; }
};

class UnaryWithProvTGSegmentItr : public TGSegmentImplItr<std::pair<Term_t,Term_t>> {
    public:
        UnaryWithProvTGSegmentItr(const size_t nodeId,
                const std::vector<std::pair<Term_t,Term_t>> &tuples) :
            TGSegmentImplItr(nodeId, tuples) {}

        Term_t get(const int colIdx) {
            return tuples[currentIdx].first;
        }

        int getNFields() const { return 1; }

        size_t getNodeId() const {
            return tuples[currentIdx].second;
        }
};

class BinaryTGSegmentItr : public TGSegmentImplItr<std::pair<Term_t,Term_t>> {
    public:
        BinaryTGSegmentItr(const size_t nodeId,
                const std::vector<std::pair<Term_t,Term_t>> &tuples) :
            TGSegmentImplItr(nodeId, tuples) {}

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
                const std::vector<BinWithProv> &tuples) :
            TGSegmentImplItr(nodeId, tuples) {}

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

#endif
