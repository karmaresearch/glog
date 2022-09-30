#include <vlog/infround/infrounditerator.h>

InfRoundIterator::InfRoundIterator(PredId_t predid,
        EDBLayer *layer,
        size_t step) :
    predid(predid), first(true), layer(layer), step(step)
{
}

Term_t InfRoundIterator::getElementAt(const uint8_t p)
{
    if (p == 0) {
        return step;
    }
    throw 10;
}
