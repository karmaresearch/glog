#ifndef _TG_CHASE_H
#define _TG_CHASE_H

#include <vlog/concepts.h>
#include <vlog/chase.h>
#include <vlog/edb.h>

struct TGChase_Node {
    Rule &rule;
};

class TGChase : public Chase {
    private:
        Program *program;
        EDBLayer &layer;

        size_t currentIteration;

#ifdef WEBINTERFACE
        PredId_t currentPredicate;
        std::string currentRule;
#endif
        std::vector<TGChase_Node> nodes;


    public:
        VLIBEXP TGChase(EDBLayer &layer, Program *program);

        VLIBEXP void run();

        Program *getProgram();

        EDBLayer &getEDBLayer();

        size_t getSizeTable(const PredId_t predid) const;

        FCIterator getTableItr(const PredId_t predid);

        FCTable *getTable(const PredId_t predid);

        size_t getCurrentIteration();

#ifdef WEBINTERFACE
        std::string getCurrentRule();

        PredId_t getCurrentPredicate();
#endif

};

#endif
