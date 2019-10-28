#include <vlog/tgchase.h>

TGChase::TGChase(EDBLayer &layer, Program *program) : layer(layer), program(program) {
}

void TGChase::run() {
    //TODO
}

Program *TGChase::getProgram() {
    return program;
}

EDBLayer &TGChase::getEDBLayer() {
    return layer;
}

size_t TGChase::getSizeTable(const PredId_t predid) const {
}

FCIterator TGChase::getTableItr(const PredId_t predid) {
}

FCTable *TGChase::getTable(const PredId_t predid) {
}

size_t TGChase::getCurrentIteration() {
}

#ifdef WEBINTERFACE

std::chrono::system_clock::time_point TGChase::getStartingTimeMs() {
}

std::string TGChase::getCurrentRule() {
}

bool TGChase::isRunning() {
}

std::vector<
std::pair<std::string, std::vector<StatsSizeIDB>>> TGChase::getSizeIDBs() {
}

std::vector<StatsRule> TGChase::getOutputNewIterations() {
}
#endif
