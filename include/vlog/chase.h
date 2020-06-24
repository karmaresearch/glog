#ifndef _CHASE_H
#define _CHASE_H

#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/fctable.h>

struct StatsSizeIDB {
    size_t step;
    int idRule;
    long derivation;
};

struct StatsRule {
    size_t step;
    long nderivations_final;
    long nderivations_unfiltered;
    long nderivations_unique;
    int idRule;

    double timems;
    double timems_first;
    double timems_merge;
    double timems_join;
    double timems_createhead;
    double timems_retain;
    std::string nbdyatoms;

    StatsRule() : idRule(-1), step(0), nderivations_final(-1),
    nderivations_unfiltered(-1), nderivations_unique(-1),
    timems(-1), timems_first(-1), timems_merge(-1),
    timems_join(-1), timems_createhead(-1), timems_retain(-1) {}
};

class Chase {
    private:
        std::string name;
        std::chrono::system_clock::time_point startTime;
        std::string profilerPath;
        bool storeStats;

#ifdef WEBINTERFACE
        long statsLastIteration;
        bool running;
#endif
        std::vector<StatsRule> statsRuleExecution;

    protected:
        virtual void saveStatistics(StatsRule &stats) {
            statsRuleExecution.push_back(stats);
        }

        void initRun();

        void stopRun();

    public:

#ifdef WEBINTERFACE
        Chase() : storeStats(true), running(false) {}
#else
        Chase() : storeStats(false) {}
#endif

        virtual Program *getProgram() = 0;

        virtual EDBLayer &getEDBLayer() = 0;

        virtual size_t getSizeTable(const PredId_t predid) const = 0;

        virtual FCIterator getTableItr(const PredId_t predid) = 0;

        virtual FCTable *getTable(const PredId_t predid) = 0;

        virtual void run() = 0;

        virtual size_t getCurrentIteration() = 0;

        void setName(const std::string &name) {
            this->name = name;
        }

        void setPathStoreStatistics(std::string profilerPath) {
            this->profilerPath = profilerPath;
            this->storeStats = true;
        }

        bool shouldStoreStats() {
            return storeStats && profilerPath != "";
        }

        const std::string &getName() const {
            return name;
        }

        std::chrono::system_clock::time_point getStartingTimeMs();

#ifdef WEBINTERFACE
        virtual std::string getCurrentRule() = 0;

        virtual PredId_t getCurrentPredicate() = 0;

        bool isRunning();

        std::vector<std::pair<std::string,
            std::vector<StatsSizeIDB>>> getSizeIDBs();

        std::vector<StatsRule> getOutputNewIterations();
#endif

};

#endif
