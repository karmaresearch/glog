#ifndef _GBSEGMENT_BINARY_H
#define _GBSEGMENT_BINARY_H

template<typename K>
bool invertedSorter(const K &a, const K &b) {
    return a.second < b.second || (a.second == b.second && a.first < b.first);
}

template<typename S, typename K, typename I, SegProvenanceType CP>
class BinaryTGSegmentImpl : public TGSegmentImpl<S,K,I,CP> {
    public:
        BinaryTGSegmentImpl(std::vector<K> &tuples,
                const size_t nodeId,
                bool isSorted=false, uint8_t sortedField = 0) :
            TGSegmentImpl<S,K,I,CP>(tuples, nodeId,
                    isSorted, sortedField) {}

        BinaryTGSegmentImpl(std::shared_ptr<std::vector<K>> tuples,
                const size_t nodeId, bool isSorted=false, uint8_t sortedField = 0)
            : TGSegmentImpl<S,K,I,CP>(tuples, nodeId, isSorted, sortedField) {}

        size_t getNColumns() const {
            return 2;
        }

        virtual std::string getName() const {
            return "TGBinarySegment";
        }

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const {
            uint8_t field = (fields.size() == 0 || fields[0] == 0) ? 0 : 1;
            std::vector<K> sortedTuples(*TGSegmentImpl<S,K,I,CP>::tuples.get());
            if (field == 0) {
                std::sort(sortedTuples.begin(), sortedTuples.end());
            } else {
                assert(field == 1);
                std::sort(sortedTuples.begin(), sortedTuples.end(), invertedSorter<K>);
            }
            return std::shared_ptr<TGSegment>(
                    new S(sortedTuples, TGSegmentImpl<S,K,I,CP>::getNodeId(),
                        true, field));
        }

        std::shared_ptr<TGSegment> swap() const {
            std::vector<K> newtuples;
            for(const K &t : *TGSegmentImpl<S,K,I,CP>::tuples.get()) {
                K newt = t;
                newt.first = t.second;
                newt.second = t.first;
                newtuples.push_back(newt);
            }
            return std::shared_ptr<TGSegment>(new S(newtuples,
                        TGSegmentImpl<S,K,I,CP>::getNodeId(), false, 0));
        }

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            if (colPos == 0) {
                for(const K &t : *TGSegmentImpl<S,K,I,CP>::tuples.get()) {
                    out.push_back(t.first);
                }
            } else {
                assert(colPos == 1);
                for(const K &t : *TGSegmentImpl<S,K,I,CP>::tuples.get()) {
                    out.push_back(t.second);
                }
            }
        }
};

class BinaryTGSegment : public BinaryTGSegmentImpl<BinaryTGSegment,
    std::pair<Term_t,Term_t>,BinaryTGSegmentItr, SEG_NOPROV>
{
    public:
        BinaryTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        BinaryTGSegment(std::shared_ptr<std::vector<std::pair<Term_t,Term_t>>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        size_t getNOffsetColumns() const {
            return 0;
        }

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            if (colPos == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(t.first);
                }
            } else {
                assert(colPos == 1);
                for(auto &t : *tuples.get()) {
                    out.push_back(t.second);
                }
            }
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.first));
                }
            } else if (colPos1 == colPos2 && colPos1 == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.first, t.first));
                }
            } else {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.second));
                }
            }
        }

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const {
            if (colIdx == 0)
                return tuples->at(rowIdx).first;
            else
                return tuples->at(rowIdx).second;
        }
};

class BinaryWithConstProvTGSegment : public BinaryTGSegmentImpl<
                                     BinaryWithConstProvTGSegment,
                                     std::pair<Term_t,Term_t>,
                                     BinaryTGSegmentItr, SEG_SAMENODE>
{
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
            } else if (colPos1 == 1 && colPos2 == 0){
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.second;
                    p.second = t.first;
                    p.node = nodeId;
                    out.push_back(p);
                }
            } else if (colPos1 == 0 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.first;
                    p.second = t.first;
                    p.node = nodeId;
                    out.push_back(p);
                }
            } else {
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.second;
                    p.second = t.second;
                    p.node = nodeId;
                    out.push_back(p);
                }
            }
        }

        size_t getNOffsetColumns() const {
            return 1;
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
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.first));
                }
            } else if (colPos1 == colPos2 && colPos1 == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.first, t.first));
                }
            } else {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.second));
                }
            }
        }

        size_t countHits(const std::vector<
                std::pair<Term_t,Term_t>> &terms,
                int column1, int column2) const {
            assert(column1 == 0 && column2 == 1);
            size_t c = 0;
            for(auto &t : terms) {
                if (std::binary_search(tuples->begin(),
                            tuples->end(), t))
                    c++;
            }
            return c;
        }

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const {
            if (colIdx == 0)
                return tuples->at(rowIdx).first;
            else
                return tuples->at(rowIdx).second;
        }
};

class BinaryWithProvTGSegment : public BinaryTGSegmentImpl<
                                BinaryWithProvTGSegment,
                                BinWithProv,
                                BinaryWithProvTGSegmentItr, SEG_DIFFNODES>
{
    private:
        static bool cmpFirstSecondTerm(const BinWithProv &a,
                const BinWithProv &b) {
            return a.first == b.first && a.second == b.second;
        }
        static bool sortNode(const BinWithProv &a,
                const BinWithProv &b) {
            return a.node < b.node || (a.node == b.node && a.first < b.first)
                || (a.node == b.node && a.first == b.first && a.second < b.second);
        }

    public:
        BinaryWithProvTGSegment(std::vector<BinWithProv> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        BinaryWithProvTGSegment(std::shared_ptr<std::vector<BinWithProv>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        size_t getNOffsetColumns() const {
            return 1;
        }

        bool isNodeConstant() const {
            return false;
        }

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            if (colPos == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.first, t.node));
                }
            } else if (colPos == 1) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.node));
                }
            } else {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            }
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.second;
                    p.second = t.first;
                    p.node = t.node;
                    out.push_back(p);
                }
            } else if (colPos1 == colPos2 && colPos1 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.first;
                    p.second = t.first;
                    p.node = t.node;
                    out.push_back(p);
                }
            } else {
                for(auto &t : *tuples.get()) {
                    BinWithProv p;
                    p.first = t.second;
                    p.second = t.second;
                    p.node = t.node;
                    out.push_back(p);
                }
            }
        }

        std::shared_ptr<const TGSegment> unique() const {
            auto t = std::vector<BinWithProv>(
                    *TGSegmentImpl<BinaryWithProvTGSegment, BinWithProv,
                    BinaryWithProvTGSegmentItr, SEG_DIFFNODES>::tuples.get());
            auto itr = std::unique(t.begin(), t.end(), cmpFirstSecondTerm);
            t.erase(itr, t.end());
            return std::shared_ptr<const TGSegment>(new BinaryWithProvTGSegment(t,
                        TGSegmentImpl<
                        BinaryWithProvTGSegment,
                        BinWithProv,
                        BinaryWithProvTGSegmentItr, SEG_DIFFNODES>::getNodeId(),
                        true, 0));
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
                        f_isSorted,
                        sortedField));
        }

        std::shared_ptr<const TGSegment> sortByProv() const {
            auto t = std::vector<BinWithProv>(
                    *TGSegmentImpl<BinaryWithProvTGSegment,
                    BinWithProv, BinaryWithProvTGSegmentItr,
                    SEG_DIFFNODES>::tuples.get());
            std::sort(t.begin(), t.end(), sortNode);
            return std::shared_ptr<const TGSegment>(new BinaryWithProvTGSegment(t,
                        TGSegmentImpl<
                        BinaryWithProvTGSegment,
                        BinWithProv,
                        BinaryWithProvTGSegmentItr, SEG_DIFFNODES>::getNodeId(),
                        false, 0));
        }

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const {
            if (colIdx == 0)
                return tuples->at(rowIdx).first;
            else
                return tuples->at(rowIdx).second;
        }
};

class BinaryWithConstNodeOffFullProvTGSegment : public BinaryTGSegmentImpl<
                                                BinaryWithConstNodeOffFullProvTGSegment,
                                                std::pair<Term_t,Term_t>, BinaryTGSegmentItr, SEG_FULLPROV>
{
    public:
        BinaryWithConstNodeOffFullProvTGSegment(std::vector<std::pair<Term_t,Term_t>> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        BinaryWithConstNodeOffFullProvTGSegment(std::shared_ptr<std::vector<std::pair<Term_t,Term_t>>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
            if (colPos == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(t.first);
                }
            } else {
                assert(colPos == 1);
                for(auto &t : *tuples.get()) {
                    out.push_back(t.second);
                }
            }
        }

        size_t getNOffsetColumns() const {
            return 2;
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.first));
                }
            } else if (colPos1 == colPos2 && colPos1 == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.first, t.first));
                }
            } else {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.second));
                }
            }
        }

        std::vector<Term_t> getRow(size_t rowIdx, bool addProv) const {
            std::vector<Term_t> out;
            out.push_back(tuples->at(rowIdx).first);
            out.push_back(tuples->at(rowIdx).second);
            if (addProv)
            {
                out.push_back(nodeId);
                out.push_back(rowIdx);
            }
            return out;
        }

        Term_t getOffsetAtRow(size_t rowIdx, size_t proofNr,
                size_t offsetColumnIdx) const {
            assert(proofNr == 0);
            assert(offsetColumnIdx == 0 && rowIdx < tuples->size());
            return rowIdx;
        }

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const {
            if (colIdx == 0)
                return tuples->at(rowIdx).first;
            else
                return tuples->at(rowIdx).second;
        }
};

class BinaryWithConstNodeFullProvTGSegment : public BinaryTGSegmentImpl<
                                             BinaryWithConstNodeFullProvTGSegment,
                                             BinWithOff,
                                             BinaryWithOffTGSegmentItr, SEG_FULLPROV>
{
    private:
        static bool cmpFirstSecondTerm(const BinWithOff &a,
                const BinWithOff &b) {
            return a.first == b.first && a.second == b.second;
        }
        static bool sortNode(const BinWithOff &a,
                const BinWithOff &b) {
            return a.off < b.off;
        }

        static bool cmp_hits(const BinWithOff &a, const BinWithOff &b) {
            return a.first < b.first || (a.first == b.first && a.second < b.second);
        }

    public:
        BinaryWithConstNodeFullProvTGSegment(std::vector<BinWithOff> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        BinaryWithConstNodeFullProvTGSegment(
                std::shared_ptr<std::vector<BinWithOff>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        size_t getNOffsetColumns() const {
            return 2;
        }

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            if (colPos == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.first, t.off));
                }
            } else if (colPos == 1) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.off));
                }
            } else {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            }
        }

        size_t countHits(const std::vector<
                std::pair<Term_t,Term_t>> &terms,
                int column1, int column2) const {
            size_t c = 0;
            for(auto &t : terms) {
                if (std::binary_search(tuples->begin(),
                            tuples->end(), BinWithOff(t.first, t.second),
                            BinaryWithConstNodeFullProvTGSegment::cmp_hits))
                    c++;
            }
            return c;
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<std::pair<Term_t, Term_t>> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.first, t.second));
                }
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    out.push_back(std::make_pair(t.second, t.first));
                }
            } else {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            }
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithOff> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithOff p;
                    p.first = t.second;
                    p.second = t.first;
                    p.off = t.off;
                    out.push_back(p);
                }
            } else if (colPos1 == colPos2 && colPos1 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithOff p;
                    p.first = t.first;
                    p.second = t.first;
                    p.off = t.off;
                    out.push_back(p);
                }
            } else {
                for(auto &t : *tuples.get()) {
                    BinWithOff p;
                    p.first = t.second;
                    p.second = t.second;
                    p.off = t.off;
                    out.push_back(p);
                }
            }
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithFullProv> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                for(auto &t : *tuples.get()) {
                    BinWithFullProv p;
                    p.first = t.first;
                    p.second = t.second;
                    p.node =  getNodeId();
                    p.prov = t.off;
                    out.push_back(p);
                }
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithFullProv p;
                    p.first = t.second;
                    p.second = t.first;
                    p.node =  getNodeId();
                    p.prov = t.off;
                    out.push_back(p);
                }
            } else if (colPos1 == colPos2 && colPos1 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithFullProv p;
                    p.first = t.first;
                    p.second = t.first;
                    p.node = getNodeId();
                    p.prov = t.off;
                    out.push_back(p);
                }
            } else {
                for(auto &t : *tuples.get()) {
                    BinWithFullProv p;
                    p.first = t.second;
                    p.second = t.second;
                    p.node = getNodeId();
                    p.prov = t.off;
                    out.push_back(p);
                }
            }
        }

        void appendTo(uint8_t colPos1,
                std::vector<UnWithFullProv> &out) const {
            assert(colPos1 == 0 || colPos1 == 1);
            auto nrows = tuples->size();
            for(size_t i = 0; i < nrows; ++i) {
                UnWithFullProv v;
                if (colPos1 == 0)
                    v.first = tuples->at(i).first;
                else
                    v.first = tuples->at(i).second;
                v.node = getNodeId();
                v.prov= tuples->at(i).off;
                out.push_back(v);
            }
        }

        std::shared_ptr<const TGSegment> unique() const {
            auto t = std::vector<BinWithOff>(
                    *TGSegmentImpl<BinaryWithConstNodeFullProvTGSegment, BinWithOff,
                    BinaryWithOffTGSegmentItr, SEG_FULLPROV>::tuples.get());
            auto itr = std::unique(t.begin(), t.end(), cmpFirstSecondTerm);
            t.erase(itr, t.end());
            return std::shared_ptr<const TGSegment>(
                    new BinaryWithConstNodeFullProvTGSegment(t,
                        TGSegmentImpl<
                        BinaryWithConstNodeFullProvTGSegment,
                        BinWithOff,
                        BinaryWithOffTGSegmentItr, SEG_FULLPROV>::getNodeId(),
                        true, 0));
        }

        /*std::shared_ptr<TGSegment> slice(size_t nodeId,
          const size_t start, const size_t end) const {
          std::vector<BinWithOff,> out(end - start);
          size_t m = 0;
          for(size_t j = start; j < end; ++j) {
          out[m] = tuples->at(j);
          m++;
          }
          return std::shared_ptr<TGSegment>(
          new BinaryWithConstNodeFullProvTGSegment(out,
          nodeId,
          f_isSorted,
          sortedField));
          }*/

        std::shared_ptr<const TGSegment> sortByProv() const {
            auto t = std::vector<BinWithOff>(
                    *TGSegmentImpl<BinaryWithConstNodeFullProvTGSegment,
                    BinWithOff, BinaryWithOffTGSegmentItr,
                    SEG_FULLPROV>::tuples.get());
            std::sort(t.begin(), t.end(), sortNode);
            return std::shared_ptr<const TGSegment>(
                    new BinaryWithConstNodeFullProvTGSegment(t,
                        TGSegmentImpl<
                        BinaryWithConstNodeFullProvTGSegment,
                        BinWithOff,
                        BinaryWithOffTGSegmentItr, SEG_FULLPROV>::getNodeId(),
                        false, 0));
        }

        std::vector<Term_t> getRow(size_t rowIdx, bool addProv) const {
            std::vector<Term_t> out;
            out.push_back(tuples->at(rowIdx).first);
            out.push_back(tuples->at(rowIdx).second);
            if (addProv)
            {
                out.push_back(nodeId);
                out.push_back(tuples->at(rowIdx).off);
            }
            return out;
        }

        Term_t getOffsetAtRow(size_t rowIdx, size_t proofNr,
                size_t offsetColumnIdx) const {
            assert(proofNr == 0);
            assert(offsetColumnIdx == 0 && rowIdx < tuples->size());
            return tuples->at(rowIdx).off;
        }

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const {
            if (colIdx == 0)
                return tuples->at(rowIdx).first;
            else
                return tuples->at(rowIdx).second;
        }
};

class BinaryWithFullProvTGSegment : public BinaryTGSegmentImpl<
                                    BinaryWithFullProvTGSegment,
                                    BinWithFullProv,
                                    BinaryWithFullProvTGSegmentItr, SEG_FULLPROV>
{
    private:
        static bool cmpFirstSecondTerm(const BinWithFullProv &a,
                const BinWithFullProv &b) {
            return a.first == b.first && a.second == b.second;
        }
        static bool sortNode(const BinWithFullProv &a,
                const BinWithFullProv &b) {
            return a.node < b.node || (a.node == b.node &&
                    a.first < b.first) || (a.node == b.node &&
                        a.first == b.first && a.second < b.second);
        }

    public:
        BinaryWithFullProvTGSegment(std::vector<BinWithFullProv> &tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        BinaryWithFullProvTGSegment(std::shared_ptr<std::vector<BinWithFullProv>> tuples,
                const size_t nodeId, bool isSorted, uint8_t sortedField) :
            BinaryTGSegmentImpl(tuples, nodeId, isSorted, sortedField) { }

        size_t getNOffsetColumns() const {
            return 2;
        }

        bool isNodeConstant() const {
            return false;
        }

        void appendTo(uint8_t colPos,
                std::vector<UnWithFullProv> &out) const {
            if (colPos == 0) {
                for(auto &t : *tuples.get()) {
                    UnWithFullProv p;
                    p.first = t.first;
                    p.node = t.node;
                    p.prov = t.prov;
                    out.push_back(p);
                }
            } else if (colPos == 1) {
                for(auto &t : *tuples.get()) {
                    UnWithFullProv p;
                    p.first = t.second;
                    p.node = t.node;
                    p.prov = t.prov;
                    out.push_back(p);
                }
            } else {
                LOG(ERRORL) << "Not possible";
                throw 10;
            }
        }

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithFullProv> &out) const {
            if (colPos1 == 0 && colPos2 == 1) {
                std::copy(tuples->begin(), tuples->end(), std::back_inserter(out));
            } else if (colPos1 == 1 && colPos2 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithFullProv p;
                    p.first = t.second;
                    p.second = t.first;
                    p.node = t.node;
                    p.prov = t.prov;
                    out.push_back(p);
                }
            } else if (colPos1 == colPos2 && colPos1 == 0) {
                for(auto &t : *tuples.get()) {
                    BinWithFullProv p;
                    p.first = t.first;
                    p.second = t.first;
                    p.node = t.node;
                    p.prov = t.prov;
                    out.push_back(p);
                }
            } else {
                for(auto &t : *tuples.get()) {
                    BinWithFullProv p;
                    p.first = t.second;
                    p.second = t.second;
                    p.node = t.node;
                    p.prov = t.prov;
                    out.push_back(p);
                }
            }
        }

        std::shared_ptr<const TGSegment> unique() const {
            auto t = std::vector<BinWithFullProv>(
                    *TGSegmentImpl<BinaryWithFullProvTGSegment, BinWithFullProv,
                    BinaryWithFullProvTGSegmentItr, SEG_FULLPROV>::tuples.get());
            auto itr = std::unique(t.begin(), t.end(), cmpFirstSecondTerm);
            t.erase(itr, t.end());
            return std::shared_ptr<const TGSegment>(
                    new BinaryWithFullProvTGSegment(t,
                        TGSegmentImpl<
                        BinaryWithFullProvTGSegment,
                        BinWithFullProv,
                        BinaryWithFullProvTGSegmentItr,
                        SEG_FULLPROV>::getNodeId(), true, 0));
        }

        std::shared_ptr<TGSegment> slice(size_t nodeId,
                const size_t start, const size_t end) const {
            std::vector<BinWithOff> out(end - start);
            size_t m = 0;
            for(size_t j = start; j < end; ++j) {
                out[m].first = tuples->at(j).first;
                out[m].second = tuples->at(j).second;
                out[m].off = tuples->at(j).prov;
                m++;
            }
            return std::shared_ptr<TGSegment>(
                    new BinaryWithConstNodeFullProvTGSegment(out,
                        nodeId,
                        f_isSorted,
                        sortedField));
        }

        std::shared_ptr<const TGSegment> sortByProv() const {
            auto t = std::vector<BinWithFullProv>(
                    *TGSegmentImpl<BinaryWithFullProvTGSegment,
                    BinWithFullProv, BinaryWithFullProvTGSegmentItr,
                    SEG_FULLPROV>::tuples.get());
            std::sort(t.begin(), t.end(), sortNode);
            return std::shared_ptr<const TGSegment>(
                    new BinaryWithFullProvTGSegment(t,
                        TGSegmentImpl<
                        BinaryWithFullProvTGSegment,
                        BinWithFullProv,
                        BinaryWithFullProvTGSegmentItr,
                        SEG_FULLPROV>::getNodeId(), false, 0));
        }

        std::vector<Term_t> getRow(size_t rowIdx, bool addProv) const {
            std::vector<Term_t> out;
            out.push_back(tuples->at(rowIdx).first);
            out.push_back(tuples->at(rowIdx).second);
            if (addProv)
            {
                out.push_back(tuples->at(rowIdx).node);
                out.push_back(tuples->at(rowIdx).prov);
            }
            return out;
        }

        Term_t getOffsetAtRow(size_t rowIdx, size_t proofNr,
                size_t offsetColumnIdx) const {
            assert(proofNr == 0);
            assert(offsetColumnIdx == 0 && rowIdx < tuples->size());
            return tuples->at(rowIdx).prov;
        }

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const {
            if (colIdx == 0)
                return tuples->at(rowIdx).first;
            else
                return tuples->at(rowIdx).second;
        }
};

#endif
