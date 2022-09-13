#include <vlog/timestamp/timestampiterator.h>

TimestampIterator::TimestampIterator(PredId_t predid, EDBLayer *layer) :
    predid(predid), first(true), layer(layer)
{
    auto t = std::chrono::system_clock::to_time_t(
            std::chrono::system_clock::now());
    auto t_str = std::ctime(&t);
    bool ok = layer->getOrAddDictNumber(t_str, strlen(t_str)-1, currentTime);
    assert(ok);
}

Term_t TimestampIterator::getElementAt(const uint8_t p)
{
    if (p == 0) {
        return currentTime;
    }
    throw 10;
}
