#include <vlog/gbchase/gbsegmentinserter.h>
#include <vlog/gbchase/gbsegment.h>

std::unique_ptr<GBSegmentInserter> GBSegmentInserter::getInserter(size_t card,
        size_t nodeColumns,
        bool delDupl) {
    if (card == 1) {
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
                    new GBSegmentInserterNAry(card));
        }
    } else {
        //singleton 
        return std::unique_ptr<GBSegmentInserter>(
                    new GBSegmentInserterNAry(card));
    }
}

void GBSegmentInserter::add(Term_t *row) {
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
