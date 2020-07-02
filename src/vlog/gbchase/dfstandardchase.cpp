#include <vlog/gbchase/dfstandardchase.h>

size_t DFStandardChase::executeRulesInStratum(
        const std::vector<size_t> &ruleIdxs,
        const size_t stratumLevel,
        const size_t stepStratum,
        size_t &step) {

    std::vector<size_t> datalogrules;
    std::vector<size_t> existentialrules;
    for(auto ruleIdx : ruleIdxs) {
        auto &rule = rules[ruleIdx];
        if (rule.isExistential()) {
            existentialrules.push_back(ruleIdx);
        } else {
            datalogrules.push_back(ruleIdx);
        }
    }

    size_t startStep = step;
    size_t derivedTuples = 0;
    //First, execute the datalog rules until saturation
    do {
        //Identify the rules that can be executed
        std::vector<std::pair<size_t,size_t>> admissibleRules;
        determineAdmissibleRules(ruleIdxs, stratumLevel, stepStratum,
                step, admissibleRules);

        std::vector<GBRuleInput> newnodes;
        for(auto p : admissibleRules) {
            auto ruleIdx = p.first;
            size_t prevstep = p.second;
            prepareRuleExecutionPlans(ruleIdx, prevstep, step, newnodes);
        }

        //Execute the rules
        auto nnodes = g.getNNodes();
        auto nodesToProcess = newnodes.size();
        LOG(INFOL) << "Nodes to process " << nodesToProcess;
        for(size_t idxNode = 0; idxNode < nodesToProcess; ++idxNode) {
            executeRule(newnodes[idxNode]);
        }

        if (nnodes != g.getNNodes()) {
            for(size_t idx = nnodes; idx < g.getNNodes(); ++idx) {
                auto nrows = g.getNodeSize(idx);
                derivedTuples += nrows;
            }
            step++;
        } else {
            break;
        }
    } while (true);

    //Now execute the existential rules one-by-one
    for(auto ruleIdx : existentialrules) {
        auto response = determineAdmissibleRule(
                ruleIdx,
                stratumLevel,
                stepStratum,
                startStep);
        const bool ruleIsAdmissible = response.first;
        if (ruleIsAdmissible) {
            size_t prevstep = response.second;
            std::vector<GBRuleInput> newnodes;
            prepareRuleExecutionPlans(ruleIdx, prevstep, step, newnodes);

            //Execute the rule associated to the node
            auto nnodes = g.getNNodes();
            auto nodesToProcess = newnodes.size();
            LOG(INFOL) << "Nodes to process " << nodesToProcess;
            for(size_t idxNode = 0; idxNode < nodesToProcess; ++idxNode) {
                executeRule(newnodes[idxNode]);
            }

            if (nnodes != g.getNNodes()) {
                for(size_t idx = nnodes; idx < g.getNNodes(); ++idx) {
                    auto nrows = g.getNodeSize(idx);
                    derivedTuples += nrows;
                }
                step++;
            }
        }
    }
    return derivedTuples;
}
