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
                s.step = itr.getCurrentIteration();
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
        if (el.step >= nextIteration && el.step < cIt) {
            out.push_back(el);
            statsLastIteration = el.step;
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
    if (shouldStoreStats()) {
        assert(profilerPath != "");
        std::ofstream ofs(profilerPath);
        ofs << "STEP,RULEID,TIME_MS,TIME_MS_FIRST,TIME_MS_MERGE,TIME_MS_JOIN,TIME_MS_HEAD,TIME_MS_RETAIN,NDER_TOTAL,NDER_UNFIL,NDER_UNIQUE,N_BDY_ATOMS" << std::endl;
        size_t idx = 0;
        for (auto &s : statsRuleExecution) {
            ofs << idx++ << ",";
            ofs << s.step << ",";
            ofs <<s.idRule << ",";
            ofs << s.timems << ",";
            ofs << s.timems_first << ",";
            ofs << s.timems_merge << ",";
            ofs << s.timems_join << ",";
            ofs << s.timems_createhead << ",";
            ofs << s.timems_retain << ",";
            ofs << s.nderivations_final << ",";
            ofs << s.nderivations_unfiltered << ",";
            ofs << s.nderivations_unique << ",";
            ofs << s.nbdyatoms << std::endl;
        }
        ofs.close();
    }
}

std::chrono::system_clock::time_point Chase::getStartingTimeMs() {
    return startTime;
}
