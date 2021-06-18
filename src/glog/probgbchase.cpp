#include <glog/probgbchase.h>
#include <iterator>

ProbGBChase::ProbGBChase(EDBLayer &layer, Program &program,
        bool optDelProofsStaticAnalysis) :
    GBChase(layer, &program, true, GBGraph::ProvenanceType::FULLPROV, false,
            false, false, true),
    optDelProofsStaticAnalysis(optDelProofsStaticAnalysis)
{
    executor = std::unique_ptr<GBRuleExecutor>(
            new GBRuleExecutor(g, layer, &program, false, true));
}

void ProbGBChase::prepareRuleExecutionPlans(
        const size_t &ruleIdx,
        const size_t prevstep,
        const size_t step,
        std::vector<GBRuleInput> &newnodes)
{
    std::vector<GBRuleInput> newNodesPerRule;
    GBChase::prepareRuleExecutionPlans(ruleIdx, prevstep, step, newNodesPerRule);
    if (newNodesPerRule.empty())
        return;

    if (optDelProofsStaticAnalysis) {
        const Rule &rule = program->getRule(ruleIdx);
        if (rule.getHeads().size() != 1) {
            LOG(WARNL) << "Static analysis on the derivation tree does not work"
                " with multiple head atoms";
            return;
        }
        const auto &h = rule.getFirstHead();
        const auto &ruleBody = rule.getBody();

        auto allVarsH = h.getAllVars();
        assert(!h.hasRepeatedVars());

        for(int64_t m = newNodesPerRule.size() - 1; m >= 0; m--) {
            auto n = newNodesPerRule[m];
            bool shouldBeRemoved = false;
            size_t bodyAtomIdx = 0;
            for(size_t i = 0; i < n.incomingEdges.size(); ++i) {
                while (bodyAtomIdx < ruleBody.size()) {
                    if (!ruleBody[bodyAtomIdx].isEDB()) {
                        break;
                    }
                    bodyAtomIdx++;
                }
                assert(bodyAtomIdx < ruleBody.size());
                if (ruleBody[bodyAtomIdx].getSharedVars(allVarsH).size() ==
                        allVarsH.size()) {
                    auto &potentialNodes = n.incomingEdges[i];
                    for(int64_t j = potentialNodes.size() - 1; j >= 0; j--) {
                        auto potentialNode = potentialNodes[j];
                        auto potentialNodeRuleIdx = g.getNodeRuleIdx(potentialNode);
                        auto potentialNodeIncEdges = g.getNodeIncomingEdges(potentialNode);
                        if (g.doesNewNodeAppearsDerivationTree(h,
                                    ruleBody[bodyAtomIdx],
                                    potentialNodeRuleIdx,
                                    potentialNodeIncEdges)) {
                            //Remove the potential node
                            potentialNodes.erase(potentialNodes.begin() + j);
                        }
                    }
                    if (potentialNodes.empty()) {
                        shouldBeRemoved = true;
                        break;
                    }
                }
                bodyAtomIdx++;
            }
            if (shouldBeRemoved) {
                newNodesPerRule.erase(newNodesPerRule.begin() + m);
            }
        }
    }

    for (auto &n : newNodesPerRule) {
        newnodes.push_back(n);
    }
}

bool ProbGBChase::executeRule(GBRuleInput &node, bool cleanDuplicates)
{
    Rule &rule = rules[node.ruleIdx];
#ifdef DEBUG
    if (rule.isEGD()) {
        throw 10; //Not supported here
    }
    LOG(INFOL) << "Executing rule " << node.ruleIdx << " "
        << rule.tostring(program, &layer);
#endif

    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
    auto outputsRule = executor->executeRule(rule, node);
    std::chrono::duration<double, std::milli> execRuntime =
        std::chrono::system_clock::now() - start;
    durationRuleExec += execRuntime;

    std::chrono::system_clock::time_point starth =
        std::chrono::system_clock::now();
    int headIdx = 0;
    bool nonempty = false;
    size_t nders = 0;
    size_t nders_un = 0;
    auto &heads = rule.getHeads();
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
                derivations = derivations->sort();
                derivations = derivations->unique();
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
