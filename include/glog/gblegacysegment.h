#ifndef LEGACY_SEGMENT_H
#define LEGACY_SEGMENT_H

#include <glog/gbsegment.h>

class TGSegmentLegacy : public TGSegment {
    private:
        const size_t nrows;
        const bool f_isSorted;
        const uint8_t sortedField;
        std::vector<std::shared_ptr<Column>> columns;
        const SegProvenanceType provenanceType;
        size_t nprovcolumns;

        bool shouldTrackProvenance() const {
            return provenanceType != SegProvenanceType::SEG_NOPROV;
        }

        bool isProvenanceAutomatic() const;

        std::shared_ptr<const TGSegment> shuffle(
                const std::vector<size_t> &idxs) const;

    public:
        TGSegmentLegacy(const std::vector<std::shared_ptr<Column>> &columns,
                size_t nrows, bool isSorted=false, uint8_t sortedField = 0,
                SegProvenanceType provenanceType = SegProvenanceType::SEG_NOPROV,
                size_t nprovcolumns = 0) :
            nrows(nrows),
            f_isSorted(isSorted),
            sortedField(sortedField),
            columns(columns),
            provenanceType(provenanceType),
            nprovcolumns(nprovcolumns)
    {
        for(auto c : columns)
            assert(c->size() == nrows);
        if (provenanceType == SEG_SAMENODE && nprovcolumns == 0)
        {
            //Add a column with a fix value
            auto newCol = std::shared_ptr<Column>(new CompressedColumn(0, nrows));
            this->columns.push_back(newCol);
            this->nprovcolumns = 1;
        }
        assert(!shouldTrackProvenance() ||  this->nprovcolumns > 0);

    }

        size_t getNRows() const {
            return nrows;
        }

        size_t getNColumns() const {
            return columns.size() - nprovcolumns;
        }

        size_t getNOffsetColumns() const {
            return nprovcolumns;
        }

        bool isEmpty() const {
            return nrows == 0;
        }

        virtual std::string getName() const {
            return "TGSegmentLegacy";
        }

        bool hasColumnarBackend() const {
            return true;
        }

        std::shared_ptr<const Column> getColumn(size_t idx) const {
            return columns[idx];
        }

        bool isSorted() const {
            return f_isSorted && sortedField == 0;
        }

        std::vector<
            std::shared_ptr<const TGSegment>> sliceByNodes(
                    size_t startNodeIdx,
                    std::vector<size_t> &provNodes) const;

        std::shared_ptr<TGSegment> slice(const size_t nodeId,
                const size_t start,
                const size_t end) const;

        std::shared_ptr<const TGSegment> slice(
                const size_t start,
                const size_t end) const;

        std::unique_ptr<TGSegmentItr> iterator(
                std::shared_ptr<const TGSegment> selfref = NULL) const;

        bool isSortedBy(std::vector<uint8_t> &fields) const;

        std::shared_ptr<const TGSegment> sort() const;

        void argsort(std::vector<size_t> &indices) const;

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const;

        std::shared_ptr<const TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const;

        std::shared_ptr<const TGSegment> sortByProv() const;

        std::shared_ptr<const TGSegment> unique() const;

        void argunique(std::vector<size_t> &idxs) const;

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const;

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const;

        void appendTo(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<std::pair<Term_t,Term_t>> &out) const;

        void appendTo(uint8_t colPos1,
                uint8_t colPos2,
                std::vector<BinWithProv> &out) const;

        void appendTo(const std::vector<int> &posFields,
                std::vector<std::vector<Term_t>> &out,
                bool withProv = false) const;

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithFullProv> &out) const;

        void appendTo(uint8_t colPos1,
                std::vector<UnWithFullProv> &out) const;

        void projectTo(const std::vector<int> &posFields,
                std::vector<std::shared_ptr<Column>> &out) const;

        std::shared_ptr<TGSegment> swap() const;

        SegProvenanceType getProvenanceType() const;

        bool isNodeConstant() const;

        size_t getNodeId() const;

        size_t countHits(const std::vector<Term_t> &terms,
                int column) const;

        size_t countHits(const std::vector<
                std::pair<Term_t,Term_t>> &terms,
                int column1, int column2) const;

        std::vector<Term_t> getRow(size_t rowIdx) const;

        Term_t getOffsetAtRow(size_t rowIdx, size_t offsetColumnIdx) const;

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const;

        ~TGSegmentLegacy();
};

#endif
