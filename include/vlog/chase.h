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
    private:
        std::string name;
        std::chrono::system_clock::time_point startTime;

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
        Chase() : running(false) {}
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
