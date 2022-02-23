#include <glog/gbsegmentinserter.h>
#include <glog/gbsegmentcomprprov.h>
#include <glog/gbsegment.h>

std::unique_ptr<GBSegmentInserter> GBSegmentInserter::getInserter(size_t card,
        size_t nodeColumns,
        bool delDupl) {
    if (card - nodeColumns == 0) {
        return std::unique_ptr<GBSegmentInserter>(
                new GBSegmentInserterNAry(card, card - nodeColumns, delDupl));
    } else if (card == 1) {
        return std::unique_ptr<GBSegmentInserter>(new GBSegmentInserterUnary(
                    delDupl));
    } else if (card == 2) {
        return std::unique_ptr<GBSegmentInserter>(
                new GBSegmentInserterBinary(delDupl));
    } else if (card > 2) {
        if (card == 4 && nodeColumns == 2 && delDupl) {
            return std::unique_ptr<GBSegmentInserter>(
                    new GBSegmentInserterBinaryWithDoubleProv(delDupl));
        } else {
            return std::unique_ptr<GBSegmentInserter>(
                    new GBSegmentInserterNAry(card, card - nodeColumns, delDupl));
        }
    } else {
        //singleton
        return std::unique_ptr<GBSegmentInserter>(
                new GBSegmentInserterNAry(card, card - nodeColumns, delDupl));
    }
}

std::shared_ptr<const TGSegment> GBSegmentInserter::compressProvNode(
        size_t nodeId, std::shared_ptr<const TGSegment> data)
{
    bool shouldCompress = false;

    if (data->getNRows() > 1 &&
            (data->getNColumns() == 1 || data->getNColumns() == 2))
    {
        auto itr = data->iterator();
        Term_t v1, v2;
        v1 = v2 = ~0ul;
        shouldCompress = true;
        while (itr->hasNext())
        {
            itr->next();
            if (v1 == ~0ul)
            {
                v1 = itr->get(0);
                v2 = itr->get(1);
            } else {
                if (itr->get(0) != v1 || itr->get(1) != v2)
                {
                    shouldCompress = false;
                    break;
                }
            }
        }
    }

    if (shouldCompress)
    {
        LOG(DEBUGL) << "Compressing the data with " << data->getNRows() << " ...";
        auto compressedData = std::shared_ptr<const TGSegment>(
                new TGSegmentProvCompr(data->getProvenanceType(), nodeId, data));
        return compressedData;
    } else
    {
        return data;
    }
}

void GBSegmentInserter::add(Term_t *row) {
    if (shouldRemoveDuplicates) {
        processedRecords++;
        if (processedRecords % 10000000 == 0)
            LOG(DEBUGL) << "Processed records: " << processedRecords;
        if (useDuplicateMap && isInMap(row)) {
            return;
        }

        if (getNRows() > checkDuplicatesAfter) {
            //Sort and remove duplicates
            auto nRemovedRecords = removeDuplicates();
            LOG(DEBUGL) << "Removed tuples " << nRemovedRecords;
            if (nRemovedRecords > 0) {
                LOG(DEBUGL) << "Diff: " << nRemovedRecords << " " <<
                    0.8 * THRESHOLD_CHECK_DUPLICATES;
                if (!useDuplicateMap &&
                        nRemovedRecords > 0.8 * THRESHOLD_CHECK_DUPLICATES) {
                    LOG(DEBUGL) << "Decided to use the map with " <<
                        getNRows() << " elements to filter out "
                        "duplicates immediately";
                    useDuplicateMap = true;
                    populateMap();
                }
                checkDuplicatesAfter = getNRows() +
                    THRESHOLD_CHECK_DUPLICATES;
            } else {
                shouldRemoveDuplicates = false;
                LOG(DEBUGL) << "No longer check!";
            }
        }
    }

    bool ok = true;
    for(auto &fn : fns) {
        if (!fn.fn(row, fn.posArgs.data())) {
            ok = false;
            break;
        }
    }
    if (ok)
        addRow(row);
}

bool GBSegmentInserterNAry::isInMap(Term_t *row) {
    if (cardCheckDuplicates == 1) {
        throw 10;
    } else if (cardCheckDuplicates == 2) {
        return s_binary.count(std::make_pair(row[0], row[1]));
    } else {
        throw 10;
        //return s.isIn(row);
    }
}

void GBSegmentInserterNAry::populateMap() {
    std::unique_ptr<Term_t[]> row = std::unique_ptr<Term_t[]>(new Term_t[cardCheckDuplicates]);
    for(size_t i = 0; i < addedRows; ++i) {
        for(size_t j = 0; j < cardCheckDuplicates; ++j) {
            row[j] = writers[j].getValue(i);
        }
        if (cardCheckDuplicates == 1) {
            throw 10;
        } else if (cardCheckDuplicates == 2) {
            s_binary.insert(std::make_pair(row[0], row[1]));
        } else {
            throw 10;
            //s.add(row.get());
        }
    }
}

size_t GBSegmentInserterNAry::removeDuplicates() {
    size_t oldSize = addedRows;
    assert(!isFinal);
    assert(columns.empty());

    //Transform the writers into columns
    std::vector<std::shared_ptr<Column>> columns;
    for(size_t j = 0; j < card; ++j) {
        columns.push_back(writers[j].getColumn());
    }

    //Sort a vector of indices
    std::vector<size_t> idxs(addedRows);
    std::iota(idxs.begin(), idxs.end(), 0);
    if (cardCheckDuplicates == 1) {
        const auto col = columns[0];
        std::sort(idxs.begin(), idxs.end(), [&col](const size_t &a, const size_t &b) { return col->getValue(a) < col->getValue(b); });
    } else if (cardCheckDuplicates == 2) {
        const auto col1 = columns[0];
        const auto col2 = columns[1];
        std::sort(idxs.begin(), idxs.end(), [&col1, &col2](const size_t &a, const size_t &b) { return col1->getValue(a) < col1->getValue(b) ||
                (col1->getValue(a) == col1->getValue(b) && col2->getValue(a) < col2->getValue(b)); });
    } else {
        throw 10; //Not supported
    }

    //Reset the writers
    addedRows = 0;
    for(size_t j = 0; j < card; ++j) {
        writers[j] = ColumnWriter();
    }
    for(size_t i = 0; i < idxs.size(); ++i) {
        auto rowIdx = idxs[i];
        if (i == 0) {
            for(size_t j = 0; j < card; ++j) {
                writers[j].add(columns[j]->getValue(rowIdx));
            }
            addedRows++;
        } else {
            auto prevRowIdx = idxs[i - 1];
            bool isDifferent = false;
            for(size_t j = 0; j < cardCheckDuplicates; ++j) {
                if (columns[j]->getValue(rowIdx) != columns[j]->getValue(prevRowIdx)) {
                    isDifferent = true;
                    break;
                }
            }
            if (isDifferent) {
                for(size_t j = 0; j < card; ++j) {
                    writers[j].add(columns[j]->getValue(rowIdx));
                }
                addedRows++;
            }
        }
    }
    return oldSize - addedRows;
}
