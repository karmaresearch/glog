#ifndef _TG_CHASE_H
#define _TG_CHASE_H

#include <vlog/gbchase/gbgraph.h>
#include <vlog/gbchase/gbchase.h>

#include <chrono>

typedef std::vector<size_t> TGChase_NodeIDs;

class TGChase_RuleIO {
    public:
        TGChase_NodeIDs out;
        std::vector<TGChase_NodeIDs> in;
};

class TGChase : public GBChase {
    private:
        std::unordered_map<size_t, size_t> mapExternalInternalIDs;

    public:
        VLIBEXP TGChase(EDBLayer &layer, Program *program);

        void run();

        void executeRule(const size_t ruleID, const TGChase_RuleIO &io);

        //Returns the number of facts that have been removed from "nodes"
        uint64_t retainNodes(const TGChase_NodeIDs &nodes);
};

#endif
