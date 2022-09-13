#include <vlog/timestamp/timestamptable.h>
#include <vlog/timestamp/timestampiterator.h>

TimestampTable::TimestampTable(PredId_t predid,
        EDBLayer *layer) : predid(predid), layer(layer) {
}

EDBIterator *TimestampTable::getIterator(const Literal &query) {
    return new TimestampIterator(predid, layer);
}

void TimestampTable::releaseIterator(EDBIterator *itr) {
    delete itr;
}
