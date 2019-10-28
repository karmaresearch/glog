#ifndef _TG_CHASE_H
#define _TG_CHASE_H

#include <vlog/concepts.h>
#include <vlog/chase.h>
#include <vlog/edb.h>

class TGChase : public Chase {
    private:
        Program *program;
        EDBLayer &layer;

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
        std::chrono::system_clock::time_point getStartingTimeMs();

        std::string getCurrentRule();

        bool isRunning();

        std::vector<
            std::pair<std::string, std::vector<StatsSizeIDB>>> getSizeIDBs();

        std::vector<StatsRule> getOutputNewIterations();
#endif

};

#endif
