#include <vlog/gbchase/gbsegmentinserter.h>
#include <vlog/gbchase/gbsegment.h>

std::unique_ptr<GBSegmentInserter> GBSegmentInserter::getInserter(size_t card) {
    if (card == 1) {
        return std::unique_ptr<GBSegmentInserter>(new GBSegmentInserterUnary());
    } else if (card == 2) {
        return std::unique_ptr<GBSegmentInserter>(
                new GBSegmentInserterBinary());
    } else if (card > 2) {
        return std::unique_ptr<GBSegmentInserter>(
                new GBSegmentInserterNAry(card));
    } else {
        throw 10; //singleton not yet supported
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
