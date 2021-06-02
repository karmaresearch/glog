#include <glog/gblegacysegment.h>
#include <glog/gbsegmentitr.h>
#include <glog/gbsegmentinserter.h>

#include <vlog/segment.h>

#include <numeric>


std::unique_ptr<TGSegmentItr> TGSegmentLegacy::iterator(
        std::shared_ptr<const TGSegment> selfref) const {
    return std::unique_ptr<TGSegmentItr>(
            new TGSegmentLegacyItr(columns, provenanceType, nprovcolumns));
}

std::vector<Term_t> TGSegmentLegacy::getRow(size_t rowIdx) const {
    std::vector<Term_t> out;
    for(auto c : columns) {
        out.push_back(c->getValue(rowIdx));
    }
    return out;
}

Term_t TGSegmentLegacy::getOffsetAtRow(size_t rowIdx,
        size_t offsetColumnIdx) const {
    assert(rowIdx < getNRows());
    assert(offsetColumnIdx < nprovcolumns - 1);
    auto c = columns[columns.size() - nprovcolumns + 1 + offsetColumnIdx];
    return c->getValue(rowIdx);
}

Term_t TGSegmentLegacy::getValueAtRow(size_t rowIdx, size_t colIdx) const {
    return columns[colIdx]->getValue(rowIdx);
}

bool TGSegmentLegacy::isProvenanceAutomatic() const {
    assert(shouldTrackProvenance());
    assert(nprovcolumns > 0);
    size_t begin = columns.size() - nprovcolumns;
    if (!columns[begin]->isConstant()) {
        return false;
    }
    for(size_t i = 1; i < nprovcolumns; i++) {
        if (columns[begin + i]->isBackedByVector()) {
            return false;
        }
    }
    return true;
}

bool TGSegmentLegacy::isSortedBy(std::vector<uint8_t> &fields) const {
    if (fields.size() != 1)
        return false;
    auto field = fields[0];
    return (f_isSorted && sortedField == field);
}

std::vector<std::shared_ptr<const TGSegment>> TGSegmentLegacy::sliceByNodes(
        size_t startNodeIdx,
        std::vector<size_t> &provNodes) const
{
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

std::shared_ptr<TGSegment> TGSegmentLegacy::sortBy(
        std::vector<uint8_t> &fields) const {
    if (f_isSorted && (fields.size() == 0 ||
                (fields.size() == 1 && fields[0] == 0))) {
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(
                    columns, nrows, true, 0, provenanceType, nprovcolumns));
    }

    if (columns.size() == 1) {
        auto column = columns[0]->sort();
        std::vector<std::shared_ptr<Column>> columns;
        columns.push_back(column);
        return std::shared_ptr<TGSegment>(
                new TGSegmentLegacy(
                    columns, column->size(), true, 0,
                    provenanceType, nprovcolumns));
    } else {
        auto nfields = columns.size();
        auto oldcols(columns);
        Segment s(nfields, oldcols);
        auto snew = s.sortBy(&fields);
        std::vector<std::shared_ptr<Column>> newColumns;
        for(int i = 0; i < nfields; ++i) {
            newColumns.push_back(snew->getColumn(i));
        }
        return std::shared_ptr<TGSegment>(
                new TGSegmentLegacy(
                    newColumns, s.getNRows(), true, fields[0],
                    provenanceType, nprovcolumns));
    }
}

std::shared_ptr<const TGSegment> TGSegmentLegacy::unique() const {
    if (!f_isSorted || sortedField != 0) {
        LOG(ERRORL) << "unique can only be called on sorted segments";
        throw 10;
    }

    std::vector<std::shared_ptr<Column>> oldcols;
    std::vector<std::shared_ptr<Column>> newcols;
    size_t nrows;

    if (columns.size() - nprovcolumns == 0) {
        if (!shouldTrackProvenance() || getNRows() == 0) {
            throw 10; //I'm not sure about what to do here, let's capture this case with an exception for now
        } else {
            //Just copy the first row
            return this->slice(0, 1);
        }
    } else if (shouldTrackProvenance() &&
            columns[columns.size() - nprovcolumns]->isConstant()) {
        size_t nfields = columns.size() - nprovcolumns;
        for(int i = 0; i < nfields; ++i)
            oldcols.push_back(columns[i]);
        for(size_t i = 1; i < nprovcolumns; ++i) {
            oldcols.push_back(columns[nfields + i]);
        }
        std::shared_ptr<Segment> s = std::shared_ptr<Segment>(
                new Segment(oldcols.size(), oldcols));
        auto retained = SegmentInserter::unique(s, nfields);
        nrows = retained->getNRows();
        for(int i = 0; i < nfields; ++i) {
            newcols.push_back(retained->getColumn(i));
        }
        newcols.push_back(std::shared_ptr<Column>(new CompressedColumn(
                        columns[nfields]->first(), nrows)));
        for(size_t i = 1; i < nprovcolumns; ++i) {
            newcols.push_back(
                    retained->getColumn(nfields + i - 1));
        }
    } else {
        size_t nfields = columns.size() - nprovcolumns;
        std::vector<std::shared_ptr<Column>> oldcols;
        for(int i = 0; i < columns.size(); ++i)
            oldcols.push_back(columns[i]);
        std::shared_ptr<Segment> s = std::shared_ptr<Segment>(
                new Segment(columns.size(), oldcols));
        auto retained = SegmentInserter::unique(s, nfields);
        for(int i = 0; i < columns.size(); ++i) {
            newcols.push_back(retained->getColumn(i));
        }
        nrows = retained->getNRows();
    }
    return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(newcols,
                nrows, true, sortedField, provenanceType, nprovcolumns));
}

std::shared_ptr<const TGSegment> TGSegmentLegacy::sort() const {
    if (!f_isSorted || sortedField != 0) {
        auto nfields = columns.size();
        std::vector<std::shared_ptr<Column>> newcols;
        if (shouldTrackProvenance()) {
            std::vector<std::shared_ptr<Column>> oldcols;
            for(int i = 0; i < columns.size() - nprovcolumns; ++i) {
                oldcols.push_back(columns[i]);
            }
            bool skipNode = false;
            if (columns[columns.size() - nprovcolumns]->isConstant()) {
                skipNode = true;
            } else {
                oldcols.push_back(columns[columns.size() - nprovcolumns]);
            }
            //Add offsets
            for(int i = 1; i < nprovcolumns; ++i) {
                oldcols.push_back(columns[columns.size() - nprovcolumns + i]);
            }
            Segment s(oldcols.size(), oldcols);
            auto news = s.sortBy(NULL);
            for(int i = 0; i < columns.size() - nprovcolumns; ++i) {
                newcols.push_back(news->getColumn(i));
            }
            if (skipNode) {
                newcols.push_back(columns[columns.size() - nprovcolumns]);
                for(size_t i = 1; i < nprovcolumns; ++i) {
                    newcols.push_back(news->
                            getColumn(columns.size() - nprovcolumns + i-1));
                }
            } else {
                for(size_t i = 0; i < nprovcolumns; ++i) {
                    newcols.push_back(news->
                            getColumn(columns.size() - nprovcolumns + i));
                }
            }
        }  else {
            auto oldcols(columns);
            Segment s(nfields, oldcols);
            auto news = s.sortBy(NULL);
            for(int i = 0; i < news->getNColumns(); ++i) {
                newcols.push_back(news->getColumn(i));
            }
        }
        return std::shared_ptr<const TGSegment>(
                new TGSegmentLegacy(newcols, nrows, true, 0,
                    provenanceType, nprovcolumns));
    } else {
        return std::shared_ptr<const TGSegment>(
                new TGSegmentLegacy(columns, nrows, true, 0,
                    provenanceType, nprovcolumns));
    }
}

void TGSegmentLegacy::argsort(std::vector<size_t> &indices) const {
    indices.resize(nrows);
    std::iota(indices.begin(), indices.end(), 0);
    if (!f_isSorted || sortedField != 0) {
        auto nfields = columns.size();
        if (shouldTrackProvenance()) {
            std::vector<std::shared_ptr<Column>> oldcols;
            for(int i = 0; i < columns.size() - nprovcolumns; ++i) {
                oldcols.push_back(columns[i]);
            }
            bool skipNode = false;
            if (columns[columns.size() - nprovcolumns]->isConstant()) {
                skipNode = true;
            } else {
                oldcols.push_back(columns[columns.size() - nprovcolumns]);
            }
            //Add offsets
            for(int i = 1; i < nprovcolumns; ++i) {
                oldcols.push_back(columns[columns.size() - nprovcolumns + i]);
            }
            Segment s(oldcols.size(), oldcols);
            s.argsort(indices);
        }  else {
            auto oldcols(columns);
            Segment s(nfields, oldcols);
            s.argsort(indices);
        }
    }
}

void TGSegmentLegacy::argunique(std::vector<size_t> &indices) const {
    size_t nrows = getNRows();
    assert(nrows > 0);
    indices.clear();
    indices.push_back(0);
    auto nfields = columns.size() - nprovcolumns;//columns.size();
    for(size_t i = 1; i < nrows; ++i) {
        bool equal = true;
        for(size_t j = 0; j < nfields; ++j) {
            if (columns[j]->getValue(i) != columns[j]->getValue(i-1)) {
                equal = false;
                break;
            }
        }
        if (!equal) {
            indices.push_back(i);
        }
    }
}

std::shared_ptr<const TGSegment> TGSegmentLegacy::shuffle(
        const std::vector<size_t> &idxs) const {
    bool changed = false;
    for(size_t i = 0; i < idxs.size(); ++i) {
        if (idxs[i] != i) {
            changed = true;
            break;
        }
    }
    if (changed) {
        std::vector<Term_t> tuples;
        auto ncolumns = columns.size();
        std::vector<std::unique_ptr<ColumnReader>> readers;
        for(size_t i = 0; i < ncolumns; ++i) {
            readers.push_back(columns[i]->getReader());
        }
        for(size_t z = 0; z < nrows; ++z) {
            for(size_t i = 0; i < ncolumns; ++i) {
                if (!readers[i]->hasNext()) {
                    throw 10;
                }
            }
            for(size_t i = 0; i < ncolumns; ++i) {
                auto v = readers[i]->next();
                tuples.push_back(v);
            }
        }
        auto inserter = GBSegmentInserter::getInserter(ncolumns, nprovcolumns,
                false);
        std::unique_ptr<Term_t[]> row = std::unique_ptr<Term_t[]>(
                new Term_t[ncolumns]);
        for(size_t z = 0; z < idxs.size(); ++z) {
            auto idx = idxs[z];
            for(size_t i = 0; i < ncolumns; ++i) {
                row[i] = tuples[idx * ncolumns + i];
            }
            inserter->add(row.get());
        }
        return inserter->getSegment(getNodeId(),
                false, 0,
                provenanceType,
                nprovcolumns);
    } else {
        return std::shared_ptr<TGSegment>(
                new TGSegmentLegacy(columns, nrows, f_isSorted, sortedField,
                    provenanceType, nprovcolumns));
    }
}

std::shared_ptr<const TGSegment> TGSegmentLegacy::sortByProv(size_t ncols,
        std::vector<size_t> &idxs,
        std::vector<size_t> &nodes) const {
    assert(shouldTrackProvenance() == true);
    assert(columns.size() > 1);
    const auto nrows = nodes.size() / ncols;
    idxs.resize(nrows);
    for(size_t i = 0; i < nrows; ++i) idxs[i] = i;
    std::stable_sort(idxs.begin(), idxs.end(), ProvSorter(nodes.data(), ncols));
    return shuffle(idxs);
}

std::shared_ptr<TGSegment> TGSegmentLegacy::slice(const size_t nodeId,
        const size_t start,
        const size_t end) const {
    std::vector<std::shared_ptr<Column>> newcols;
    auto length = end - start;
    int ncols = shouldTrackProvenance() ?
        columns.size() - nprovcolumns : columns.size();
    for(int i = 0; i < ncols; ++i) {
        if (start > 0 || end < nrows) {
            auto c = columns[i]->slice(start, end);
            newcols.push_back(c);
        } else {
            newcols.push_back(columns[i]);
        }
    }
    if (shouldTrackProvenance()) {
        assert(nprovcolumns > 0);
        newcols.push_back(std::shared_ptr<Column>(
                    new CompressedColumn(nodeId, length)));
        for(size_t i = 1; i < nprovcolumns; ++i) {
            auto c = columns[columns.size() - nprovcolumns + i]->slice(start, end);
            newcols.push_back(c);
        }
    }
    return std::shared_ptr<TGSegment>(new TGSegmentLegacy(newcols, length,
                f_isSorted, sortedField, provenanceType, nprovcolumns));
}

std::shared_ptr<const TGSegment> TGSegmentLegacy::slice(
        const size_t start,
        const size_t end) const {
    std::vector<std::shared_ptr<Column>> newcols;
    auto length = end - start;
    for(int i = 0; i < columns.size(); ++i) {
        if (start > 0 || end < nrows) {
            auto c = columns[i]->slice(start, end);
            newcols.push_back(c);
        } else {
            newcols.push_back(columns[i]);
        }
    }
    return std::shared_ptr<TGSegment>(new TGSegmentLegacy(newcols, length,
                f_isSorted, sortedField, provenanceType, nprovcolumns));
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
    assert(shouldTrackProvenance());
    assert(provenanceType != SEG_FULLPROV);
    auto &c1 = columns[colPos];
    auto itr1 = c1->getReader();
    auto &c2 = columns[columns.size() - 1]; //Last column is the provenance
    auto itrProv = c2->getReader();
    auto nrows = c1->size();
    for(size_t i = 0; i < nrows; ++i) {
        if (!itrProv->hasNext())
            throw 10;
        out.push_back(std::make_pair(itr1->next(), itrProv->next()));
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
    assert(shouldTrackProvenance());
    assert(provenanceType != SEG_FULLPROV);
    assert(columns.size() > colPos1);
    assert(columns.size() > colPos2);
    auto &c1 = columns[colPos1];
    auto &c2 = columns[colPos2];
    auto &c3 = columns.back();
    auto itr1 = c1->getReader();
    auto itr2 = c2->getReader();
    auto itrProv = c3->getReader();
    auto nrows = c1->size();
    for(size_t i = 0; i < nrows; ++i) {
        BinWithProv v;
        v.first = itr1->next();
        v.second = itr2->next();
        if (!itrProv->hasNext()) {
            throw 10;
        }
        v.node = itrProv->next();
        out.push_back(v);
    }
}

void TGSegmentLegacy::appendTo(uint8_t colPos1, uint8_t colPos2,
        std::vector<BinWithFullProv> &out) const {
    assert(shouldTrackProvenance());
    assert(provenanceType == SEG_FULLPROV);
    assert(columns.size() > colPos1);
    assert(columns.size() > colPos2);
    auto &c1 = columns[colPos1];
    auto &c2 = columns[colPos2];
    auto &node = columns[columns.size() - nprovcolumns];
    auto &prov = columns.back();
    auto itr1 = c1->getReader();
    auto itr2 = c2->getReader();
    auto itrNode = node->getReader();
    auto itrProv = prov->getReader();
    auto nrows = c1->size();
    for(size_t i = 0; i < nrows; ++i) {
        BinWithFullProv v;
        v.first = itr1->next();
        v.second = itr2->next();
        if (!itrNode->hasNext()) {
            throw 10;
        }
        if (!itrProv->hasNext()) {
            throw 10;
        }
        v.node = itrNode->next();
        v.prov= itrProv->next();
        out.push_back(v);
    }
}

void TGSegmentLegacy::appendTo(uint8_t colPos1,
        std::vector<UnWithFullProv> &out) const {
    assert(shouldTrackProvenance());
    assert(provenanceType == SEG_FULLPROV);
    assert(columns.size() > colPos1);
    auto &c1 = columns[colPos1];
    auto &node = columns[columns.size() - nprovcolumns];
    auto &prov = columns.back();
    auto itr1 = c1->getReader();
    auto itrNode = node->getReader();
    auto itrProv = prov->getReader();
    auto nrows = c1->size();
    for(size_t i = 0; i < nrows; ++i) {
        UnWithFullProv v;
        v.first = itr1->next();
        if (!itrNode->hasNext()) {
            throw 10;
        }
        if (!itrProv->hasNext()) {
            throw 10;
        }
        v.node = itrNode->next();
        v.prov= itrProv->next();
        out.push_back(v);
    }
}

void TGSegmentLegacy::appendTo(const std::vector<int> &posFields,
        std::vector<std::vector<Term_t>> &out,
        bool withProv) const {
    if (provenanceType == SEG_FULLPROV) {
        size_t copyNProvColumns = nprovcolumns;
        if (!withProv) {
            copyNProvColumns = 0;
            assert(out.size() == posFields.size());
        } else {
            assert(out.size() >= posFields.size() + copyNProvColumns);
        }
        assert(posFields.size() > 0 || copyNProvColumns > 0); //Otherwise, we won't be able to copy anything
        size_t i = 0;
        for(auto pos : posFields) {
            auto c = columns[pos];
            auto itr = c->getReader();
            while (itr->hasNext()) {
                out[i].push_back(itr->next());
            }
            i++;
        }
        for(size_t j = 0; j < copyNProvColumns; ++j) {
            auto c = columns[columns.size() - nprovcolumns + j];
            auto itr = c->getReader();
            while (itr->hasNext()) {
                out[i].push_back(itr->next());
            }
            i++;
        }
    } else {
        if (!withProv) {
            assert(posFields.size() == out.size());
        } else {
            assert(posFields.size() == out.size() - 1);
        }
        size_t i = 0;
        for(auto pos : posFields) {
            auto c = columns[pos];
            auto itr = c->getReader();
            while (itr->hasNext()) {
                out[i].push_back(itr->next());
            }
            i++;
        }
        if (shouldTrackProvenance()) { //Add the nodeID
            auto &o = out.back();
            auto itr = columns.back()->getReader();
            if (posFields.size() > 0) {
                while (itr->hasNext()) {
                    o.push_back(itr->next());
                }
            } else {
                if (itr->hasNext()) {
                    o.push_back(itr->next());
                }
            }
        }
    }
}

std::shared_ptr<TGSegment> TGSegmentLegacy::swap() const {
    if (columns.size() - nprovcolumns != 2) {
        LOG(ERRORL) << "Not supposed to be invoked on non binary predicates";
        throw 10;
    }
    if (shouldTrackProvenance()) {
        std::vector<std::shared_ptr<Column>> c;
        c.push_back(columns[1]);
        c.push_back(columns[0]);
        for(size_t i = 0; i < nprovcolumns; ++i) {
            c.push_back(columns[columns.size() - nprovcolumns + i]);
        }
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(
                    c, nrows, false, 0, provenanceType, nprovcolumns));
    } else {
        std::vector<std::shared_ptr<Column>> c;
        c.push_back(columns[1]);
        c.push_back(columns[0]);
        return std::shared_ptr<TGSegment>(new TGSegmentLegacy(
                    c, nrows, false, 0, provenanceType, nprovcolumns));
    }
}

void TGSegmentLegacy::projectTo(const std::vector<int> &fields,
        std::vector<std::shared_ptr<Column>> &out) const {
    for (auto &field : fields)
        out.push_back(columns[field]);
    if (shouldTrackProvenance()) {
        for(size_t i = 0; i < nprovcolumns; ++i) {
            out.push_back(columns[columns.size() - nprovcolumns + i]);
        }
    }
}

SegProvenanceType TGSegmentLegacy::getProvenanceType() const {
    assert(provenanceType != SEG_SAMENODE || columns[columns.size() - nprovcolumns]->isConstant());
    assert(provenanceType != SEG_DIFFNODES || !columns[columns.size() - nprovcolumns]->isConstant());
    return provenanceType;
}

size_t TGSegmentLegacy::getNodeId() const {
    if (shouldTrackProvenance()) {
        assert(columns.size() > 0);
        if (!columns[columns.size() - nprovcolumns]->isEmpty() &&
                columns[columns.size() - nprovcolumns]->isConstant()) {
            return columns[columns.size() - nprovcolumns]->first();
        } else {
            return ~0ul;
        }
    } else {
        return ~0ul;
    }
}

std::shared_ptr<const TGSegment> TGSegmentLegacy::sortByProv() const {
    assert(shouldTrackProvenance() == true);
    assert(columns.size() > 1);
    //if (!isProvenanceAutomatic()) {
    std::vector<uint8_t> sortedFields;
    //Node ID
    sortedFields.push_back(columns.size() - nprovcolumns);
    //Other fields
    for(int i = 0; i < columns.size() - nprovcolumns; ++i) {
        if (i != columns.size() - nprovcolumns)
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
    return std::shared_ptr<TGSegment>(
            new TGSegmentLegacy(newcols, nrows, true, 0,
                provenanceType, nprovcolumns));
    //} else {
    //    return std::shared_ptr<TGSegment>(
    //            new TGSegmentLegacy(columns, nrows, f_isSorted,
    //                sortedField, provenanceType, nprovcolumns));
    //}


    /*if (!shouldTrackProvenance()) {
      LOG(ERRORL) << "This method should not be called if the segment"
      " does not support the provenance";
      throw 10;
      }
      auto nfields = columns.size();
      auto oldcols(columns);
      Segment s(nfields, oldcols);
      std::vector<uint8_t> fields;
      fields.push_back(nfields - 1);
      auto snew = s.sortBy(&fields);
      std::vector<std::shared_ptr<Column>> newColumns;
      for(int i = 0; i < nfields; ++i) {
      newColumns.push_back(snew->getColumn(i));
      }
      return std::shared_ptr<TGSegment>(
      new TGSegmentLegacy(newColumns, s.getNRows(), false, fields[0],
      provenanceType));*/
}

size_t TGSegmentLegacy::countHits(const std::vector<Term_t> &terms,
        int column) const {
    return columns[column]->countHits(terms);
}

bool TGSegmentLegacy::isNodeConstant() const {
    assert(nprovcolumns > 0);
    return columns[columns.size() - nprovcolumns]->isConstant();
}

size_t TGSegmentLegacy::countHits(const std::vector<std::pair<
        Term_t,Term_t>> &terms,
        int column1, int column2) const {
    auto c1 = getColumn(column1);
    auto c2 = getColumn(column2);
    if (c1->isEDB() && c2->isEDB()) {
        EDBColumn *ec1 = (EDBColumn*)c1.get();
        EDBColumn *ec2 = (EDBColumn*)c2.get();
        const Literal &l1 = ec1->getLiteral();
        const Literal &l2 = ec2->getLiteral();
        std::vector<Substitution> subs;
        if (!l1.sameVarSequenceAs(l2) ||
                l1.subsumes(subs, l1, l2) == -1) {
            //The columns come from different literals. This is not supported
            throw 10;
        }
        EDBLayer &layer = ec1->getEDBLayer();
        auto v = layer.checkNewIn(terms, l1, ec1->posColumnInLiteral(),
                ec2->posColumnInLiteral());
        return terms.size() - v.size();
    } else {
        assert(column1 == 0 && column2 == 1); //I must ensure the columns are sorted
        size_t count = 0;
        auto itr1 = c1->getReader();
        auto itr2 = c2->getReader();
        size_t countTerms = 0;
        while (itr1->hasNext() && itr2->hasNext()) {
            if (countTerms == terms.size())
                break;
            auto v1 = itr1->next();
            auto v2 = itr2->next();
            if (v1 < terms[countTerms].first ||
                    (v1 == terms[countTerms].first && v2 < terms[countTerms].second)) {
                //Go to the next one
            } else {
                while (countTerms < terms.size()) {
                    if (v1 == terms[countTerms].first && v2 == terms[countTerms].second) {
                        count++;
                    } else if (v1 < terms[countTerms].first || (v1 == terms[countTerms].first && v2 < terms[countTerms].second)) {
                        break;
                    }
                    countTerms++;
                }
            }
        }
        return count;
    }
}

TGSegmentLegacy::~TGSegmentLegacy() {
}
