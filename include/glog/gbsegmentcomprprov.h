#ifndef _GBSEGMENT_PROVCOMPR_H
#define _GBSEGMENT_PROVCOMPR_H

#include <glog/gbsegment.h>

class TGSegmentProvCompr : public TGSegment {
    private:
        SegProvenanceType provenanceType;
        size_t nodeId;
        uint32_t arity;
        std::vector<Term_t> data;

        uint32_t arity_provenance;
        std::vector<std::pair<size_t, size_t>> provenance_startend;
        std::vector<size_t> provenance_offsets;

    public:

        TGSegmentProvCompr(
                const SegProvenanceType provenanceType,
                size_t nodeId,
                std::shared_ptr<const TGSegment> segment);

        size_t getNRows() const;

        size_t getNColumns() const;

        size_t getNOffsetColumns() const;

        bool isEmpty() const;

        std::unique_ptr<TGSegmentItr> iterator(
                std::shared_ptr<const TGSegment> selfref = NULL) const;

        bool isSortedBy(std::vector<uint8_t> &fields) const;

        bool isSorted() const;

        std::shared_ptr<const TGSegment> sort() const;

        void argsort(std::vector<size_t> &idxs) const;

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const;

        std::shared_ptr<const TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const;

        std::shared_ptr<const TGSegment> sortByProv() const;

        std::shared_ptr<const TGSegment> unique() const;

        void argunique(std::vector<size_t> &idxs) const;

        void projectTo(const std::vector<int> &posFields,
                std::vector<std::shared_ptr<Column>> &out) const;

        Term_t getOffsetAtRow(size_t rowIdx,
                size_t proofNr, size_t offsetColumnIdx) const;

        size_t getNProofsAtRow(size_t rowIdx) const;

        std::vector<Term_t> getRow(size_t rowIdx, bool addProv) const;

        //0 = no provenance, 1 = all tuples come from the same node
        //2 = tuples from different nodes
        SegProvenanceType getProvenanceType() const;

        std::shared_ptr<const TGSegment> slice(
                const size_t start,
                const size_t end) const;

        virtual size_t getNodeId() const;

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithFullProv> &out) const;

        void appendTo(uint8_t colPos, std::vector<UnWithFullProv> &out) const;

        Term_t getValueAtRow(size_t rowIdx, size_t colIdx) const;
};
#endif
