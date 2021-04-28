#ifndef _SEGMENT_CONCATENATED_H
#define _SEGMENT_CONCATENATED_H

#include <vlog/gbchase/gbsegment.h>
#include <vlog/gbchase/gbgraph.h>

class CompositeTGSegment : public TGSegment {
    private:
        const size_t nodeId;
        const GBGraph &g;
        std::vector<size_t> nodes;
        std::vector<int> copyVarPos;
        bool f_isSorted;
        size_t sortedField;
        const bool trackProvenance;

        const bool sortBeforeAccess;
        const bool removeDuplBeforeAccess;

        std::shared_ptr<const TGSegment> merge() const;

    public:
        CompositeTGSegment(size_t nodeId, const GBGraph &g,
                const std::vector<size_t> &nodes,
                const std::vector<int> &copyVarPos,
                bool isSorted=false, uint8_t sortedField = 0,
                bool trackProvenance = false,
                bool sortBeforeAccess = false,
                bool removeDuplBeforeAccess = false) :
            nodeId(nodeId),
            g(g),
            nodes(nodes),
            copyVarPos(copyVarPos),
            f_isSorted(isSorted),
            sortedField(sortedField),
            trackProvenance(trackProvenance),
            sortBeforeAccess(sortBeforeAccess),
            removeDuplBeforeAccess(removeDuplBeforeAccess)
    {
        assert(trackProvenance);
#ifdef DEBUG
                //There should not be repeated variables in copyVarPos, I'm
                //not sure it works if there are
                std::set<size_t> s;
                for(auto c : copyVarPos)
                    s.insert(c);
                if (s.size() != copyVarPos.size()) {
                    LOG(ERRORL) << "Case not implemented";
                    throw 10;
                }
#endif
    }

        CompositeTGSegment(const GBGraph &g,
                const std::vector<size_t> &nodes,
                const std::vector<int> &copyVarPos,
                bool isSorted=false, uint8_t sortedField = 0,
                bool trackProvenance = false,
                bool sortBeforeAccess = false,
                bool removeDuplBeforeAccess = false) :
            CompositeTGSegment(~0ul, g, nodes, copyVarPos, isSorted, sortedField,
                    trackProvenance, sortBeforeAccess, removeDuplBeforeAccess) {
            }

        std::string getName() const {
            return "CompositeTGSegment";
        }

        size_t getNColumns() const {
            return copyVarPos.size();
        }

        bool isSorted() const {
            return f_isSorted && sortedField == 0;
        }

        bool isEmpty() const {
            return false;
        }

        SegProvenanceType getProvenanceType() const {
            if (trackProvenance) {
                if (nodes.size() == 1)
                    return SEG_SAMENODE;
                else
                    return SEG_DIFFNODES;
            } else
                return SEG_NOPROV;
        }

        size_t getNodeId() const {
            return nodeId;
        }

        size_t getNRows() const;

        std::shared_ptr<TGSegment> swap() const;

        std::unique_ptr<TGSegmentItr> iterator(
                std::shared_ptr<const TGSegment> selfref = NULL) const;

        bool isSortedBy(std::vector<uint8_t> &fields) const;

        std::shared_ptr<const TGSegment> sort() const;

        std::shared_ptr<TGSegment> sortBy(std::vector<uint8_t> &fields) const;

        std::shared_ptr<TGSegment> sortByProv(size_t ncols,
                std::vector<size_t> &idxs,
                std::vector<size_t> &nodes) const;

        std::shared_ptr<const TGSegment> sortByProv() const;

        std::shared_ptr<const TGSegment> unique() const;

        void projectTo(const std::vector<int> &posFields,
                std::vector<std::shared_ptr<Column>> &out) const;

        std::vector<std::shared_ptr<const TGSegment>> sliceByNodes(
                size_t startNodeIdx,
                std::vector<size_t> &provNodes) const;

        void appendTo(uint8_t colPos, std::vector<Term_t> &out) const;

        void appendTo(uint8_t colPos,
                std::vector<std::pair<Term_t, Term_t>> &out) const;

        void appendTo(uint8_t colPos1, uint8_t colPos2,
                std::vector<BinWithProv> &out) const;
};

#endif
