#ifndef _ELASTIC_ITR_H
#define _ELASTIC_ITR_H

#include <vlog/column.h>
#include <vlog/edbtable.h>
#include <vlog/edbiterator.h>
#include <vlog/segment.h>

class ElasticIterator : public EDBIterator {
    private:
        PredId_t predid;
        uint64_t t1, t2;
        int64_t idx;

    public:
        ElasticIterator(PredId_t predid, uint64_t t1, uint64_t t2) :
            predid(predid), t1(t1), t2(t2), idx(-1) {
            }

        bool hasNext() {
            return idx == -1;
        }

        void next() {
            idx++;
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
            idx = -1;
        }
};

#endif
