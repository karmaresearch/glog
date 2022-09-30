#include <vlog/infround/infroundtable.h>
#include <vlog/infround/infrounditerator.h>
#include <vlog/edb.h>

InfRoundTable::InfRoundTable(PredId_t predid,
        EDBLayer *layer) : predid(predid), layer(layer)
{
}

void InfRoundTable::setContext(GBGraph *g, size_t step)
{
    std::string sStep = std::to_string(step);
    layer->getOrAddDictNumber(sStep.c_str(), sStep.size(), currentStep, true);
}

EDBIterator *InfRoundTable::getIterator(const Literal &query)
{
    return new InfRoundIterator(predid, layer, currentStep);
}

void InfRoundTable::releaseIterator(EDBIterator *itr)
{
    delete itr;
}
