#include <vlog/tgchase.h>

TGChase::TGChase(EDBLayer &layer, Program *program) : layer(layer), program(program) {
}

Program *TGChase::getProgram() {
    return program;
}

EDBLayer &TGChase::getEDBLayer() {
    return layer;
}

void TGChase::run() {
    initRun();
    size_t nnodes = 0;
    do {
        nnodes = nodes.size();
        //TODO: Execute the rules
    } while (nnodes != nodes.size());

    stopRun();
}

size_t TGChase::getSizeTable(const PredId_t predid) const {
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

FCIterator TGChase::getTableItr(const PredId_t predid) {
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

FCTable *TGChase::getTable(const PredId_t predid) {
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

size_t TGChase::getCurrentIteration() {
    return currentIteration;
}

#ifdef WEBINTERFACE
std::string TGChase::getCurrentRule() {
    return currentRule;
}

PredId_t TGChase::getCurrentPredicate() {
    return currentPredicate;
}
#endif
