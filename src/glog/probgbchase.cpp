#include <glog/probgbchase.h>

ProbGBChase::ProbGBChase(EDBLayer &layer, Program &program) :
    GBChase(layer, &program, true, GBGraph::ProvenanceType::FULLPROV, false,
            false, false, true)
{
    executor = std::unique_ptr<GBRuleExecutor>(
            new GBRuleExecutor(g, layer, &program, false, true));
}

bool ProbGBChase::executeRule(GBRuleInput &node, bool cleanDuplicates)
{
    Rule &rule = rules[node.ruleIdx];
#ifdef DEBUG
    if (rule.isEGD()) {
        throw 10; //Not supported here
    }
    LOG(DEBUGL) << "Executing rule " << node.ruleIdx << " "
        << rule.tostring(program, &layer);
#endif

    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();

    size_t nders = 0;
    size_t nders_un = 0;
    auto &heads = rule.getHeads();
    int headIdx = 0;
    auto outputsRule = executor->executeRule(rule, node);
    bool nonempty = false;
    std::chrono::duration<double, std::milli> execRuntime =
        std::chrono::system_clock::now() - start;
    durationRuleExec += execRuntime;

    std::chrono::system_clock::time_point starth =
        std::chrono::system_clock::now();
    for (GBRuleOutput &outputRule : outputsRule) {
        auto currentPredicate = heads[headIdx++].getPredicate().getId();
        auto derivations = outputRule.segment;
        nders_un += derivations->getNRows();
        triggers += derivations->getNRows();
        auto derivationNodes = outputRule.nodes;
        nonempty |= !(derivations == NULL || derivations->isEmpty());

        if (nonempty) {
            nders += derivations->getNRows();
            assert(heads.size() == 1);
            bool filterProvenance = !heads.back().isMagic();
            if (filterProvenance) {
                g.addNodesProv(currentPredicate, node.ruleIdx,
                        node.step, derivations, derivationNodes, true);
            } else {
                auto retainedTuples = g.retain(currentPredicate, derivations,
                        derivationNodes);
                nonempty = !(retainedTuples == NULL || retainedTuples->isEmpty());
                if (nonempty) {
                    g.addNodesProv(currentPredicate, node.ruleIdx,
                            node.step, retainedTuples, derivationNodes, false);
                }
            }
        }
    }

    if (shouldStoreStats()) {
        std::chrono::duration<double, std::milli> retainRuntime =
            std::chrono::system_clock::now() - starth;
        std::chrono::duration<double, std::milli> totalRuntime =
            std::chrono::system_clock::now() - start;

        StatsRule stats;
        stats.step = node.step;
        stats.idRule = node.ruleIdx;
        stats.nderivations_final = nders;
        stats.nderivations_unfiltered = nders_un;
        stats.timems = totalRuntime.count();
        stats.timems_first = executor->getDuration(DurationType::DUR_FIRST).count();
        stats.timems_merge = executor->getDuration(DurationType::DUR_MERGE).count();
        stats.timems_join = executor->getDuration(DurationType::DUR_JOIN).count();
        stats.timems_createhead = executor->getDuration(DurationType::DUR_HEAD).count();
        stats.timems_retain = retainRuntime.count();
        stats.nbdyatoms = executor->getStat(StatType::N_BDY_ATOMS);
        saveStatistics(stats);
    }

    return nonempty;
}
