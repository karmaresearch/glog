#ifndef _CLIQUE_ITERATOR_H
#define _CLIQUE_ITERATOR_H

#include <vlog/column.h>
#include <vlog/edbtable.h>
#include <vlog/edbiterator.h>
#include <vlog/segment.h>

class CliqueIterator : public EDBIterator {
    private:
        PredId_t predid;
        const bool empty;
        const std::map<size_t, std::vector<Term_t>> &components;

        std::map<Term_t, size_t>::iterator begitr;
        std::map<Term_t, size_t>::iterator itr;
        std::map<Term_t, size_t>::iterator enditr;

        std::vector<Term_t>::const_iterator startitrComp;
        std::vector<Term_t>::const_iterator enditrComp;

        bool hasAdvanced;
        bool advance();
        Term_t t1, t2;

    public:
        CliqueIterator(PredId_t predid,
                const std::map<size_t, std::vector<Term_t>> &components,
                std::map<Term_t, size_t>::iterator startitr,
                std::map<Term_t, size_t>::iterator enditr);

        bool hasNext();

        void next();

        Term_t getElementAt(const uint8_t p);

        void skipDuplicatedFirstColumn();

        void clear();

        PredId_t getPredicateID() {
            return predid;
        }
};

#endif
