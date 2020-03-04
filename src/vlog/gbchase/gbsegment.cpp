#include <vlog/gbchase/gbsegment.h>
#include <vlog/gbchase/gbsegmentitr.h>

#include <vlog/segment.h>

std::unique_ptr<TGSegmentItr> TGSegmentLegacy::iterator() const {
    return std::unique_ptr<TGSegmentItr>(new TGSegmentLegacyItr(columns, trackProvenance));
}

bool TGSegmentLegacy::isSortedBy(std::vector<uint8_t> &fields) const {
    if (fields.size() != 1)
        return false;
    auto field = fields[0];
    return (isSorted && sortedField == field);
}

std::shared_ptr<TGSegment> TGSegmentLegacy::sortBy(std::vector<uint8_t> &fields) const {
    if (isSorted && (fields.size() == 0 || (fields.size() == 1 && fields[0] == 0))) {
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(columns, nrows, true, 0, trackProvenance));
    }
    if (columns.size() == 1) {
        auto column = columns[0]->sort();
        std::vector<std::shared_ptr<Column>> columns;
        columns.push_back(column);
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(columns, column->size(), true, 0, trackProvenance));
    } else {
        auto nfields = columns.size();
        auto oldcols(columns);
        Segment s(nfields, oldcols);
        auto snew = s.sortBy(&fields);
        std::vector<std::shared_ptr<Column>> newColumns;
        for(int i = 0; i < nfields; ++i) {
            newColumns.push_back(snew->getColumn(i));
        }
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(newColumns, s.getNRows(), true, fields[0], trackProvenance));
    }
}

std::shared_ptr<TGSegment> TGSegmentLegacy::unique() const {
    if (!isSorted || sortedField != 0) {
        LOG(ERRORL) << "unique can only be called on sorted segments";
        throw 10;
    }
    auto nfields = columns.size();
    auto nfieldsTocheck = trackProvenance ? nfields - 1 : nfields;
    auto oldcols(columns);
    std::shared_ptr<Segment> s = std::shared_ptr<Segment>(new Segment(nfields, oldcols));
    auto retained = SegmentInserter::unique(s, nfieldsTocheck);
    std::vector<std::shared_ptr<Column>> newcols;
    for(int i = 0; i < retained->getNColumns(); ++i) {
        newcols.push_back(retained->getColumn(i));
    }
    return std::shared_ptr<TGSegment>(new TGSegmentLegacy(newcols, retained->getNRows(), true, sortedField, trackProvenance));
}

std::shared_ptr<TGSegment> TGSegmentLegacy::sort() const {
    if (!isSorted || sortedField != 0) {
        auto nfields = columns.size();
        auto oldcols(columns);
        Segment s(nfields, oldcols);
        auto news = s.sortBy(NULL);
        std::vector<std::shared_ptr<Column>> newcols;
        for(int i = 0; i < news->getNColumns(); ++i) {
            newcols.push_back(news->getColumn(i));
        }
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(newcols, nrows, true, 0, trackProvenance));
    } else {
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(columns, nrows, true, 0, trackProvenance));
    }
}

std::shared_ptr<TGSegment> TGSegmentLegacy::sortByProv(size_t ncols,
        std::vector<size_t> &idxs,
        std::vector<size_t> &nodes) const {
    assert(trackProvenance == true);
    assert(columns.size() > 1);
    if (!columns.back()->isConstant()) {
        std::vector<uint8_t> sortedFields;
        sortedFields.push_back(columns.size() - 1);
        for(int i = 0; i < columns.size() - 1; ++i) {
            sortedFields.push_back(i);
        }
        auto nfields = columns.size();
        auto oldcols(columns);
        Segment s(nfields, oldcols);
        auto news = s.sortBy(&sortedFields);
        std::vector<std::shared_ptr<Column>> newcols;
        for(int i = 0; i < news->getNColumns(); ++i) {
            newcols.push_back(news->getColumn(i));
        }
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(newcols, nrows, true, 0, trackProvenance));
    } else {
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(columns, nrows, isSorted, sortedField, trackProvenance));
    }
}

std::shared_ptr<TGSegment> TGSegmentLegacy::slice(const size_t nodeId,
        const size_t start,
        const size_t end) const {
    std::vector<std::shared_ptr<Column>> newcols;
    auto length = end - start;
    int ncols = trackProvenance ? columns.size() - 1 : columns.size();
    for(int i = 0; i < ncols; ++i) {
        if (start > 0 || end < nrows - 1) {
            auto &v = columns[i]->getVectorRef();
            std::vector<Term_t> slicedColumn(length);
            std::copy(v.begin() + start, v.begin() + end, slicedColumn.begin());
            std::shared_ptr<Column> c(new InmemoryColumn(slicedColumn, true));
            newcols.push_back(c);
        } else {
            newcols.push_back(columns[i]);
        }
    }
    if (trackProvenance) {
        CompressedColumnBlock b(nodeId, 0, length);
        std::vector<CompressedColumnBlock> blocks;
        blocks.push_back(b);
        newcols.push_back(std::shared_ptr<Column>(
                    new CompressedColumn(blocks, length)));
    }
    return std::shared_ptr<TGSegment>(new TGSegmentLegacy(newcols, length, isSorted, sortedField, trackProvenance));

}

void TGSegmentLegacy::appendTo(uint8_t colPos, std::vector<Term_t> &out) const {
    assert(colPos < columns.size());
    auto &col = columns[colPos];
    auto itr = col->getReader();
    while (itr->hasNext()) {
        out.push_back(itr->next());
    }
}

void TGSegmentLegacy::appendTo(uint8_t colPos,
        std::vector<std::pair<Term_t, Term_t>> &out) const {
    assert(trackProvenance);
    auto &c1 = columns[colPos];
    auto &v1 = c1->getVectorRef();
    auto &c2 = columns[columns.size() - 1]; //Last column is the provenance
    auto itrProv = c2->getReader();
    auto nrows = v1.size();
    for(size_t i = 0; i < nrows; ++i) {
        if (!itrProv->hasNext())
            throw 10;
        out.push_back(std::make_pair(v1[i], itrProv->next()));
    }
}

void TGSegmentLegacy::appendTo(uint8_t colPos1,
        uint8_t colPos2,
        std::vector<std::pair<Term_t,Term_t>> &out) const {
    auto &c1 = columns[colPos1];
    auto &c2 = columns[colPos2];
    auto itr1 = c1->getReader();
    auto itr2 = c2->getReader();
    auto nrows = c1->size();
    for(size_t i = 0; i < nrows; ++i) {
        out.push_back(std::make_pair(itr1->next(), itr2->next()));
    }
}

void TGSegmentLegacy::appendTo(uint8_t colPos1,
        uint8_t colPos2,
        std::vector<BinWithProv> &out) const {
    assert(trackProvenance);
    auto &c1 = columns[colPos1];
    auto &c2 = columns[colPos2];
    auto &c3 = columns[2];
    auto &v1 = c1->getVectorRef();
    auto &v2 = c2->getVectorRef();
    auto itrProv = c3->getReader();
    auto nrows = v1.size();
    for(size_t i = 0; i < nrows; ++i) {
        BinWithProv v;
        v.first = v1[i];
        v.second = v2[i];
        if (!itrProv->hasNext()) {
            throw 10;
        }
        v.node = itrProv->next();
        out.push_back(v);
    }
}

std::shared_ptr<TGSegment> TGSegmentLegacy::swap() const {
    if (trackProvenance) {
        if (columns.size() != 3) {
            LOG(ERRORL) << "Not supposed to be invoked on non binary predicates";
            throw 10;
        }
        std::vector<std::shared_ptr<Column>> c;
        c.push_back(columns[1]);
        c.push_back(columns[0]);
        c.push_back(columns[2]);
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(
                    c, nrows, false, 0, trackProvenance));
    } else {
        if (columns.size() != 2) {
            LOG(ERRORL) << "Not supposed to be invoked on non binary predicates";
            throw 10;
        }
        std::vector<std::shared_ptr<Column>> c;
        c.push_back(columns[1]);
        c.push_back(columns[0]);
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(
                    c, nrows, false, 0, trackProvenance));
    }
}

void TGSegmentLegacy::projectTo(const std::vector<int> &fields,
        std::vector<std::shared_ptr<Column>> &out) const {
    for (auto &field : fields)
        out.push_back(columns[field]);
    if (trackProvenance) {
        out.push_back(columns.back());
    }
}

int TGSegmentLegacy::getProvenanceType() const {
    if (trackProvenance) {
        assert(columns.size() > 0);
        if (columns.back()->isConstant()) {
            return 1;
        } else {
            return 2;
        }
    } else {
        return 0;
    }
}

size_t TGSegmentLegacy::getNodeId() const {
    if (trackProvenance) {
        assert(columns.size() > 0);
        assert(columns.back()->size() > 0);
        if (columns.back()->isConstant()) {
            return columns.back()->first();
        } else {
            return ~0ul;
        }
    } else {
        return ~0ul;
    }
}
