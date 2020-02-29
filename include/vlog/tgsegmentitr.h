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
        virtual bool hasNext() = 0;

        virtual void next() = 0;

        virtual Term_t get(const int colIdx) = 0;

        virtual size_t getNodeId() const = 0;

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

        int getNFields() const {
            return trackProvenance ? columns.size() - 1 : columns.size();
        }
};

/*class TGSegmentLegacyDirectItr : public TGSegmentDirectItr {
    private:
        const bool trackProvenance;
        const std::vector<std::shared_ptr<Column>> &columns;

        int64_t nrows;
        std::vector<const Term_t*> vectors;
        std::vector<Term_t> constValues;
        int64_t currentRowIdx;
        int64_t m_currentRowIdx;

    public:
        TGSegmentLegacyDirectItr(const std::vector<std::shared_ptr<Column>> &columns,
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
                    auto c = columns[i];
                    if (c->isConstant()) {
                        vectors.push_back(NULL);
                        constValues.push_back(c->first());
                    } else {
                        vectors.push_back(c->getVectorRef().data());
                        constValues.push_back(0);
                    }
                }
            }

        bool hasNext() {
            return currentRowIdx < (nrows - 1);
        }

        void next() {
            currentRowIdx++;
        }

        void mark() {
            m_currentRowIdx = currentRowIdx;
        }

        void reset() {
            currentRowIdx = m_currentRowIdx;
        }

        Term_t get(const int colIdx) {
            if (!vectors[colIdx]) {
                return constValues[colIdx];
            } else {
                return vectors[colIdx][currentRowIdx];
            }
        }

        size_t getNodeId() const {
            return columns.back()->getValue(currentRowIdx);
        }

        int getNFields() const {
            return trackProvenance ? columns.size() - 1 : columns.size();
        }
};*/

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
