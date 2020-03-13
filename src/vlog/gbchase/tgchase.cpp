#include <vlog/gbchase/tgchase.h>

TGChase::TGChase(EDBLayer &layer, Program *program) : GBChase(layer, program) {
}

void TGChase::run() {
    LOG(ERRORL) << "Run should not invoked with TGChase";
    throw 10;
}

void TGChase::executeRule(const size_t ruleID, const TGChase_RuleIO &io) {
    //TODO
}

uint64_t TGChase::retainNodes(const TGChase_NodeIDs &nodes) {
    std::vector<size_t> internalIDs;
    for(auto node : nodes) {
        assert(mapExternalInternalIDs.count(node));
        internalIDs.push_back(mapExternalInternalIDs[node]);
    }
    return g.removeDuplicatesFromNodes(internalIDs);
}
