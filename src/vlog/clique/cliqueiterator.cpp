#include <vlog/clique/cliqueiterator.h>

CliqueIterator::CliqueIterator(PredId_t predid,
        const std::map<size_t, std::vector<Term_t>> &components,
        std::map<Term_t, size_t>::iterator startitr,
        std::map<Term_t, size_t>::iterator enditr) :
    predid(predid), components(components), startitr(startitr), enditr(enditr) {

        auto compId = startitr->second;
        assert(components.count(compId));
        const auto &comp = components.find(compId);
        const auto &vector = comp->second;
        assert(vector.size() > 0);
        startitrComp = vector.cbegin();
        enditrComp = vector.cend();
        hasAdvanced = true;
        t1 = t2 = ~0ul;
    }

bool CliqueIterator::advance() {
    if (startitrComp != enditrComp) {
        assert(startitrComp < enditrComp);
        startitrComp++;
    }
    if (startitrComp == enditrComp) {
        //Move to the next term
        startitr++;
        if (startitr != enditr) {
            auto compId = startitr->second;
            assert(components.count(compId));
            const auto &comp = components.find(compId);
            startitrComp = comp->second.cbegin();
            enditrComp = comp->second.cend();
        } else {
            return false;
        }
    } else {
    }
    return true;
}

bool CliqueIterator::hasNext() {
    if (hasAdvanced) {
        return true;
    } else {
        hasAdvanced = advance();
        return hasAdvanced;
    }
}

void CliqueIterator::next() {
    if (!hasAdvanced) {
        hasAdvanced = advance();
        if (!hasAdvanced) {
            t1 = t2 = ~0ul;
            return;
        }
    }
    if (!hasAdvanced || startitrComp == enditrComp)
        throw 10;
    t1 = startitr->first;
    t2 = *startitrComp;
    hasAdvanced = false;
}

Term_t CliqueIterator::getElementAt(const uint8_t p) {
    if (p == 0) {
        return t1;
    } else {
        assert(p == 1);
        return t2;
    }
}

void CliqueIterator::skipDuplicatedFirstColumn() {
    LOG(ERRORL) << "Not supported";
    throw 10;
}

void CliqueIterator::clear() {
    LOG(ERRORL) << "Not supported";
    throw 10;
}
