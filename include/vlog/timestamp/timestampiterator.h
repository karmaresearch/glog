#ifndef _TIMESTAMP_ITR_H
#define _TIMESTAMP_ITR_H

#include <vlog/edbiterator.h>
#include <vlog/edb.h>
#include <chrono>

class TimestampIterator : public EDBIterator {
    private:
        PredId_t predid;
        EDBLayer *layer;
        bool first;
        Term_t currentTime;

    public:
        TimestampIterator(PredId_t predid, EDBLayer *layer);
                
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
