#ifndef _STRING_ITR_H
#define _STRING_ITR_H

#include <vlog/column.h>
#include <vlog/edbtable.h>
#include <vlog/edbiterator.h>
#include <vlog/segment.h>

class StringIterator : public EDBIterator {
    private:
        PredId_t predid;
        uint64_t t1, t2;
        bool found;
        bool processed;

    public:
        StringIterator(PredId_t predid, uint64_t t1, uint64_t t2,
                bool found) :
            predid(predid), t1(t1), t2(t2), found(found), processed(false) {
            }

        StringIterator(PredId_t predid, uint64_t t1,
                bool found) :
            predid(predid), t1(t1), found(found), processed(false) {
            }


        bool hasNext() {
            return found && !processed;
        }

        void next() {
            processed = true;
        }

        Term_t getElementAt(const uint8_t p) {
            if (p == 0) {
                return t1;
            } else {
                return t2;
            }
        }

        PredId_t getPredicateID() {
            return predid;
        }

        void skipDuplicatedFirstColumn() {
        }

        void clear() {
            processed = false;
        }
};

#endif
