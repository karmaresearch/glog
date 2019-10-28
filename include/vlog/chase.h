#ifndef _CHASE_H
#define _CHASE_H

#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/fctable.h>

struct StatsSizeIDB {
    size_t iteration;
    int idRule;
    long derivation;
};

struct StatsRule {
    size_t iteration;
    long derivation;
    int idRule;
    long timems;
    long totaltimems;
    StatsRule() : idRule(-1) {}
};

class Chase {

    public:
        virtual Program *getProgram() = 0;
        virtual EDBLayer &getEDBLayer() = 0;

        virtual size_t getSizeTable(const PredId_t predid) const = 0;
        virtual FCIterator getTableItr(const PredId_t predid) = 0;
        virtual FCTable *getTable(const PredId_t predid) = 0;

        virtual void run() = 0;

        virtual size_t getCurrentIteration() = 0;

#ifdef WEBINTERFACE

        virtual std::chrono::system_clock::time_point getStartingTimeMs() = 0;

        virtual std::string getCurrentRule() = 0;

        virtual bool isRunning() = 0;

        virtual std::vector<
            std::pair<std::string, std::vector<StatsSizeIDB>>> getSizeIDBs() = 0;

        virtual std::vector<StatsRule> getOutputNewIterations() = 0;
#endif

};

#endif
