#include <glog/gbsegmentcomprprov.h>
#include <glog/gbsegmentitr.h>

TGSegmentProvCompr::TGSegmentProvCompr(
        SegProvenanceType provenanceType,
        size_t nodeId,
        std::shared_ptr<const TGSegment> segment)
{
    this->provenanceType = provenanceType;
    this->nodeId = nodeId;
    assert(nodeId != ~0ul);
    this->arity = segment->getNColumns();
    assert(arity > 0);
    //The segment must be sorted
    assert(segment->isSorted());
    this->arity_provenance = segment->getNOffsetColumns() - 1; //The first one
    //is the node
    auto itr = segment->iterator();
    std::vector<Term_t> vs(this->arity);
    for(auto &v : vs)
        v = ~0ul;
    while (itr->hasNext())
    {
        itr->next();
        bool nw = false;
        for (size_t i = 0; i < this->arity; ++i)
            if (vs[i] != itr->get(i))
            {
                nw = true;
                break;
            }
        if (nw)
        {
            for (size_t i = 0; i < this->arity; ++i)
            {
                data.push_back(itr->get(i));
            }
            provenance_startend.push_back(std::make_pair(
                        provenance_startend.size(),
                        provenance_startend.size()));
            for (size_t i = 0; i < this->arity; ++i)
                vs[i] = itr->get(i);
        }
        //Copy the provenance
        assert(itr->getNProofs() == 1);
        for(size_t i = 0; i < this->arity_provenance; ++i)
        {
            provenance_offsets.push_back(itr->getProvenanceOffset(0, i));
        }
        provenance_startend.back().second += this->arity_provenance;
    }
}

size_t TGSegmentProvCompr::getNRows() const
{
    return data.size() / arity;
}

size_t TGSegmentProvCompr::getNColumns() const
{
    return arity;
}

size_t TGSegmentProvCompr::getNOffsetColumns() const
{
    return arity_provenance + 1;
}

bool TGSegmentProvCompr::isEmpty() const
{
    return data.empty();
}

std::vector<Term_t> TGSegmentProvCompr::getRow(size_t rowIdx, bool addProv) const
{
    assert(addProv == false);
    std::vector<Term_t> out;
    size_t startIdx = rowIdx * arity;
    for(int i = startIdx; i < startIdx + arity; ++i)
    {
        out.push_back(data[i]);
    }
    return out;
}

size_t TGSegmentProvCompr::getNProofsAtRow(size_t rowIdx) const
{
    return (provenance_startend[rowIdx].second -
            provenance_startend[rowIdx].first) / arity_provenance;
}

Term_t TGSegmentProvCompr::getOffsetAtRow(size_t rowIdx,
        size_t proofNr, size_t offsetColumnIdx) const
{
    auto range = provenance_startend[rowIdx];
    auto idx = range.first + proofNr * arity_provenance + offsetColumnIdx;
    assert(idx < range.second);
    return provenance_offsets[idx];
}

std::unique_ptr<TGSegmentItr> TGSegmentProvCompr::iterator(
        std::shared_ptr<const TGSegment> selfref) const
{
    return std::unique_ptr<TGSegmentItr>(
            new TGSegmentCompProvItr(provenanceType, nodeId, arity,
                data));
}

bool TGSegmentProvCompr::isSortedBy(std::vector<uint8_t> &fields) const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

bool TGSegmentProvCompr::isSorted() const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

std::shared_ptr<const TGSegment> TGSegmentProvCompr::sort() const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

void TGSegmentProvCompr::argsort(std::vector<size_t> &idxs) const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

std::shared_ptr<TGSegment> TGSegmentProvCompr::sortBy(std::vector<uint8_t> &fields) const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

std::shared_ptr<const TGSegment> TGSegmentProvCompr::sortByProv(size_t ncols,
        std::vector<size_t> &idxs,
        std::vector<size_t> &nodes) const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

std::shared_ptr<const TGSegment> TGSegmentProvCompr::sortByProv() const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

std::shared_ptr<const TGSegment> TGSegmentProvCompr::unique() const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

void TGSegmentProvCompr::argunique(std::vector<size_t> &idxs) const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

void TGSegmentProvCompr::projectTo(const std::vector<int> &posFields,
        std::vector<std::shared_ptr<Column>> &out) const
{
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

//0 = no provenance, 1 = all tuples come from the same node
//2 = tuples from different nodes
SegProvenanceType TGSegmentProvCompr::getProvenanceType() const
{
    return provenanceType;
}

std::shared_ptr<const TGSegment> TGSegmentProvCompr::slice(
        const size_t start,
        const size_t end) const
{
    throw 10;
}

size_t TGSegmentProvCompr::getNodeId() const
{
    return nodeId;
}

void TGSegmentProvCompr::appendTo(uint8_t colPos1, uint8_t colPos2,
        std::vector<BinWithFullProv> &out) const
{
    auto itr = iterator();
    size_t counter = 0;
    while (itr->hasNext())
    {
        itr->next();
        BinWithFullProv v;
        v.first = itr->get(colPos1);
        v.second = itr->get(colPos2);
        v.node = nodeId;
        v.prov= counter++;
        out.push_back(v);
    }
}

void TGSegmentProvCompr::appendTo(uint8_t colPos,
        std::vector<UnWithFullProv> &out) const
{
    auto itr = iterator();
    size_t counter = 0;
    while (itr->hasNext())
    {
        itr->next();
        UnWithFullProv v;
        v.first = itr->get(colPos);
        v.node = nodeId;
        v.prov= counter++;
        out.push_back(v);
    }
}
