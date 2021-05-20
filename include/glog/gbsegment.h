#ifndef _TG_SEGMENT_H
#define _TG_SEGMENT_H

#include <vlog/concepts.h>
#include <vlog/column.h>
#include <vlog/segment.h>

#include <glog/gbsegmentitr.h>

#include <vector>
#include <map>
#include <numeric>

class TGSegment {
    private:

    public:
        virtual size_t getNRows() const = 0;

        virtual size_t getNColumns() const = 0;

        virtual bool isEmpty() const = 0;

        virtual std::unique_ptr<TGSegmentItr> iterator(
                std::shared_ptr<const TGSegment> selfref = NULL) const = 0;

        virtual bool isSortedBy(std::vector<uint8_t> &fields) const = 0;

        virtual bool isSorted() const = 0;

        virtual std::shared_ptr<const TGSegment> sort() const = 0;

        virtual void argsort(std::vector<size_t> &idxs) const = 0;

        virtual std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const = 0;

        virtual std::shared_ptr<const TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const = 0;

        virtual std::shared_ptr<const TGSegment> sortByProv() const = 0;

        virtual std::shared_ptr<const TGSegment> unique() const = 0;

        virtual void argunique(std::vector<size_t> &idxs) const = 0;

        virtual void projectTo(const std::vector<int> &posFields,
                std::vector<std::shared_ptr<Column>> &out) const = 0;

        //0 = no provenance, 1 = all tuples come from the same node
        //2 = tuples from different nodes
        virtual SegProvenanceType getProvenanceType() const = 0;

        virtual size_t getNodeId() const = 0;

        virtual size_t countHits(const std::vector<Term_t> &terms,
                int column) const {
            LOG(ERRORL) << "(Hits) Not implemented";
            throw 10;
        }

        virtual size_t countHits(const std::vector<
                std::pair<Term_t,Term_t>> &terms,
                int column1, int column2) const {
            LOG(ERRORL) << "(Hits) Not implemented";
            throw 10;
        }

        virtual std::shared_ptr<TGSegment> slice(const size_t nodeId,
                const size_t start,
                const size_t end) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual std::vector<
            std::shared_ptr<const TGSegment>> sliceByNodes(
                    size_t startNodeIdx,
                    std::vector<size_t> &provNodes) const {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            }

        virtual std::shared_ptr<const TGSegment> shuffle(
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

        virtual void appendTo(uint8_t colPos,
                std::vector<UnWithFullProv> &out) const {
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

        virtual void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithFullProv> &out) const {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }

        virtual void appendTo(const std::vector<int> &posFields,
                std::vector<std::vector<Term_t>> &out,
                bool withProv = false) const {
            if (posFields.size() == 0) {
                assert(out.size() == 1);
                out[0].push_back(getNodeId());
            } else {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            }
        }

        virtual void projectTo(uint8_t colPos, std::vector<Term_t> &out) const {
            out.clear();
            appendTo(colPos, out);
        }

        virtual void projectTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            out.clear();
            appendTo(colPos, out);
        }

        virtual void projectTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            out.clear();
            appendTo(colPos1, colPos2, out);
        }

        virtual void projectTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            out.clear();
            appendTo(colPos1, colPos2, out);
        }

        virtual bool hasColumnarBackend() const {
            return false;
        }

        virtual size_t getNOffsetColumns() const {
            return 0;
        }

        virtual bool isNodeConstant() const {
            return true;
        }

        virtual ~TGSegment() {}
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
            if (tuples[ba + i] < tuples[bb + i])
                return true;
            else if (tuples[ba + i] > tuples[bb + i]) {
                return false;
            }
        }
        return false;
    }
};

//S is the name of the class, K is the undelying container,
//I is the iterator, CP is a flag that indicates whether the container
//supports provenance (0=no provenance, 1=constant provenance,
//2=variable provenance
template<typename S, typename K, typename I, SegProvenanceType CP>
class TGSegmentImpl : public TGSegment {
    protected:
        std::shared_ptr<std::vector<K>> tuples;
        const size_t nodeId;
        const bool f_isSorted;
        const uint8_t sortedField;

    public:
        TGSegmentImpl(std::vector<K> &t, const size_t nodeId, bool isSorted=false,
                uint8_t sortedField=0) :
            nodeId(nodeId), f_isSorted(isSorted), sortedField(sortedField) {
                tuples = std::shared_ptr<std::vector<K>>(new std::vector<K>());
                tuples->swap(t);
            }

        TGSegmentImpl(std::shared_ptr<std::vector<K>> t, const size_t nodeId, bool isSorted=false,
                uint8_t sortedField=0) :
            tuples(t), nodeId(nodeId), f_isSorted(isSorted), sortedField(sortedField) {
            }

        size_t getNRows() const {
            return tuples->size();
        }

        const std::vector<K> &getTuples() const {
            return *tuples.get();
        }

        size_t getNodeId() const {
            return nodeId;
        }

        bool isEmpty() const {
            return tuples->empty();
        }

        SegProvenanceType getProvenanceType() const {
            return CP;
        }

        virtual std::shared_ptr<const TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const {
            const auto nrows = nodes.size() / ncols;
            idxs.resize(nrows);
            for(size_t i = 0; i < nrows; ++i) idxs[i] = i;
            std::stable_sort(idxs.begin(), idxs.end(), ProvSorter(nodes.data(), ncols));
            return shuffle(idxs);
        }

        virtual std::shared_ptr<const TGSegment> sortByProv() const {
            LOG(ERRORL) << "Not supported";
            throw 10;
        }

        virtual std::shared_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) const {
            if (start == 0 && end == getNRows()) {
                return std::shared_ptr<TGSegment>(new S(tuples,
                            nodeId,
                            TGSegmentImpl<S,K,I,CP>::f_isSorted,
                            TGSegmentImpl<S,K,I,CP>::sortedField));
            } else {
                std::vector<K> out(end - start);
                std::copy(TGSegmentImpl<S,K,I,CP>::tuples->begin() + start,
                        TGSegmentImpl<S,K,I,CP>::tuples->begin() + end, out.begin());
                return std::shared_ptr<TGSegment>(new S(out,
                            nodeId,
                            TGSegmentImpl<S,K,I,CP>::f_isSorted,
                            TGSegmentImpl<S,K,I,CP>::sortedField));
            }
        }

        std::shared_ptr<const TGSegment> shuffle(const std::vector<size_t> &idxs) const {
            std::vector<K> out;
            for(auto idx : idxs) out.push_back(TGSegmentImpl<S,K,I,CP>::tuples->at(idx));
            return std::shared_ptr<const TGSegment>(new S(out, TGSegmentImpl<S,K,I,CP>::getNodeId(), false, 0));
        }

        std::unique_ptr<TGSegmentItr> iterator(
                std::shared_ptr<const TGSegment> selfref = NULL) const {
            return std::unique_ptr<TGSegmentItr>(new I(
                        TGSegmentImpl<S,K,I,CP>::getNodeId(),
                        *TGSegmentImpl<S,K,I,CP>::tuples.get(),
                        selfref));
        }

        bool isSortedBy(std::vector<uint8_t> &fields) const {
            if (fields.size() != 1)
                return false;
            auto field = fields[0];
            assert(field == 0 || field == 1);
            return TGSegmentImpl<S,K,I,CP>::f_isSorted && field == TGSegmentImpl<S,K,I,CP>::sortedField;
        }

        bool isSorted() const {
            return f_isSorted && sortedField == 0;
        }

        std::shared_ptr<const TGSegment> sort() const {
            auto t = std::vector<K>(*TGSegmentImpl<S,K,I,CP>::tuples.get());
            std::sort(t.begin(), t.end());
            return std::shared_ptr<const TGSegment>(new S(t, TGSegmentImpl<S,K,I,CP>::getNodeId(), true, 0));
        }

        void argsort(std::vector<size_t> &indices) const {
            const auto t = std::vector<K>(*TGSegmentImpl<S,K,I,CP>::tuples.get());
            indices.resize(t.size());
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(),
                    [&t](size_t left, size_t right) -> bool {
                    return t[left] < t[right];
                    });
        }

        virtual std::shared_ptr<const TGSegment> unique() const {
            auto t = std::vector<K>(*TGSegmentImpl<S,K,I,CP>::tuples.get());
            auto itr = std::unique(t.begin(), t.end());
            t.erase(itr, t.end());
            return std::shared_ptr<const TGSegment>(new S(t, TGSegmentImpl<S,K,I,CP>::getNodeId(), true, 0));
        }

        virtual void argunique(std::vector<size_t> &idxs) const {
            const auto t = std::vector<K>(*TGSegmentImpl<S,K,I,CP>::tuples.get());
            assert(!t.empty());
            idxs.clear();
            idxs.push_back(0);
            for(size_t i = 1; i < t.size(); ++i) {
                if (t[i] == t[i-1]) {
                } else {
                    idxs.push_back(i);
                }
            }
        }

        virtual void projectTo(const std::vector<int> &posFields,
                std::vector<std::shared_ptr<Column>> &out) const {
            assert(posFields.size() > 2);
            out.clear();
            auto itr = iterator();
            const int nfields = CP == 2 ? posFields.size() + 1 : posFields.size();
            SegmentInserter ins(nfields);
            size_t count = 0;
            if (CP == 2) {
                while (itr->hasNext()) {
                    itr->next();
                    for(int i = 0; i < posFields.size(); ++i) {
                        ins.addAt(i, itr->get(posFields[i]));
                    }
                    ins.addAt(posFields.size(), itr->getNodeId());
                    count++;
                }
            } else {
                while (itr->hasNext()) {
                    itr->next();
                    for(int i = 0; i < posFields.size(); ++i) {
                        ins.addAt(i, itr->get(posFields[i]));
                    }
                    count++;
                }
            }
            auto seg = ins.getSegment();
            for(int i = 0; i < posFields.size(); ++i) {
                out.push_back(seg->getColumn(i));
            }
            if (CP == 1) {
                out.push_back(std::shared_ptr<Column>(
                            new CompressedColumn(getNodeId(), count)));
            } else if (CP == 2) {
                out.push_back(seg->getColumn(posFields.size()));
            }
        }

        virtual std::vector<std::shared_ptr<const TGSegment>> sliceByNodes(
                size_t startNodeIdx,
                std::vector<size_t> &provNodes) const {
            std::vector<std::shared_ptr<const TGSegment>> out;
            size_t startidx = 0;
            auto itr = iterator();
            size_t i = 0;
            size_t currentNode = ~0ul;
            while (itr->hasNext()) {
                itr->next();
                bool hasChanged = i == 0 || currentNode != itr->getNodeId();
                if (hasChanged) {
                    if (startidx < i) {
                        //Create a new node
                        provNodes.push_back(currentNode);
                        auto dataToAdd = slice(
                                startNodeIdx++, startidx, i);
                        out.push_back(dataToAdd);
                    }
                    startidx = i;
                    currentNode = itr->getNodeId();
                }
                i++;
            }
            //Copy the last segment
            if (startidx < i) {
                auto dataToAdd = slice(startNodeIdx++, startidx, i);
                provNodes.push_back(currentNode);
                out.push_back(dataToAdd);
            }
            return out;
        }


        virtual ~TGSegmentImpl() {}
};

#include <glog/gbsegment_unary.h>

#include <glog/gbsegment_binary.h>


#endif
