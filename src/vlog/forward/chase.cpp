#include <vlog/chase.h>
#include <vlog/ruleexecdetails.h>

#ifdef WEBINTERFACE
std::vector<std::pair<std::string, std::vector<StatsSizeIDB>>> Chase::getSizeIDBs() {
    std::vector<std::pair<std::string, std::vector<StatsSizeIDB>>> out;
    Program *program = getProgram();
    for (PredId_t i = 0; i < program->getNPredicates(); ++i) {
        if (i != getCurrentPredicate() && program->isPredicateIDB(i)) {
            FCIterator itr = getTableItr(i);
            std::vector<StatsSizeIDB> stats;
            while (!itr.isEmpty()) {
                std::shared_ptr<const FCInternalTable> t = itr.getCurrentTable();
                StatsSizeIDB s;
                s.iteration = itr.getCurrentIteration();
                s.idRule = itr.getRule()->ruleid;
                s.derivation = t->getNRows();
                stats.push_back(s);
                itr.moveNextCount();
            }
            if (stats.size() > 0) {
                out.push_back(std::make_pair(program->getPredicateName(i), stats));
            }
        }
    }
    return out;
}

std::vector<StatsRule> Chase::getOutputNewIterations() {
    std::vector<StatsRule> out;
    size_t cIt = getCurrentIteration();
    int nextIteration = statsLastIteration + 1;
    size_t sizeVector = statsRuleExecution.size();
    for (int i = 0; i < sizeVector; ++i) {
        StatsRule el = statsRuleExecution[i];
        if (el.iteration >= nextIteration && el.iteration < cIt) {
            out.push_back(el);
            statsLastIteration = el.iteration;
        }
    }
    return out;
}

bool Chase::isRunning() {
    return running;
}
#endif

void Chase::initRun() {
#ifdef WEBINTERFACE
    statsLastIteration = -1;
    running = true;
#endif
    startTime = std::chrono::system_clock::now();
    statsRuleExecution.clear();
}

void Chase::stopRun() {
#ifdef WEBINTERFACE
    running = false;
#endif
}

std::chrono::system_clock::time_point Chase::getStartingTimeMs() {
    return startTime;
}
