#ifndef _INFROUND_ITR_H
#define _INFROUND_ITR_H

#include <vlog/edbiterator.h>
#include <vlog/edb.h>
#include <chrono>

class InfRoundIterator : public EDBIterator {
    private:
        PredId_t predid;
        EDBLayer *layer;
        bool first;
        size_t step;

    public:
        InfRoundIterator(PredId_t predid, EDBLayer *layer, size_t step);
                
        bool hasNext() {
            return first;
        }

        void next() {
            first = false;
        }

        Term_t getElementAt(const uint8_t p);

        PredId_t getPredicateID() {
            return predid;
        }

        void skipDuplicatedFirstColumn() {
        }

        void clear() {
            first = true;
        }
};

#endif
