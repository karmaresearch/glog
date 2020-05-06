#include <vlog/seminaiver_trigger.h>
#include <vlog/gbchase/tgchase.h>
#include <vlog/trigger/tgpath.h>

TGChase::TGChase(EDBLayer &layer, Program *program) : GBChase(layer, program) {
}

void TGChase::executeRule(const size_t ruleIdx, TGChase_RuleIO &io) {
    if (trackProvenance) {
        size_t nOutputNodes = io.inNodes[0].size();
        for(size_t i = 1; i < io.inNodes.size(); ++i) {
            nOutputNodes = nOutputNodes * io.inNodes[i].size();
        }
        if (nOutputNodes != io.outNodes.size()) {
            LOG(ERRORL) << "The number of output nodes does not match the "
                "combination of nodes in input";
            throw 10;
        }
    } else {
        if (io.outNodes.size() != 1) {
            LOG(ERRORL) << "The configuration did not specify the output "
                "nodes";
            throw 10;
        }
    }

    GBRuleInput input;
    input.ruleIdx = ruleIdx;
    input.step = 1;
    for (auto &nodes : io.inNodes)  {
        std::vector<size_t> nodeIDs;
        for (auto node : nodes) {
            assert(mapExternalInternalIDs.count(node));
            nodeIDs.push_back(mapExternalInternalIDs[node]);
        }
        input.incomingEdges.push_back(nodeIDs);
    }

    size_t nNodes = g.getNNodes();
    bool newNodes = GBChase::executeRule(input, false);
    io.outSizes.clear();
    if (newNodes) {
        if (trackProvenance) {
            LOG(ERRORL) << "Not implemented (yet)";
            throw 10;
        } else {
            //There should be only one node added
            size_t nNodesAfterExecution = g.getNNodes();
            if (nNodesAfterExecution != nNodes + 1) {
                LOG(ERRORL) << "There should have been only one atom added";
                throw 10;
            }
            std::string nodeOut = io.outNodes[0];
            size_t sizeNodeOut = g.getNodeSize(nNodes);
            mapExternalInternalIDs[nodeOut] = nNodes;
            io.outSizes.push_back(sizeNodeOut);
        }
    } else {
        for(auto &outNode : io.outNodes) {
            io.outSizes.push_back(0);
        }
    }
}

/*uint64_t TGChase::retainNodes(const TGChase_NodeIDs &nodes) {
  std::vector<size_t> internalIDs;
  for(auto node : nodes) {
  assert(mapExternalInternalIDs.count(node));
  internalIDs.push_back(mapExternalInternalIDs[node]);
  }
  return g.removeDuplicatesFromNodes(internalIDs);
  }*/

void TGChaseStatic::run() {
    initRun();
    LOG(DEBUGL) << "First load the trigger_path file";
    TGPaths paths(tgfile);
    LOG(DEBUGL) << "There are " << paths.getNPaths() << " paths to execute";

    //Read the file one path at the time and execute the rules
    std::set<PredId_t> createdPredicates;
    std::unordered_map<std::string, size_t> iterations;
    std::vector<StatIteration> costRules;
    size_t nTriggers = 0;

    for(size_t i = 0; i < paths.getNPaths(); ++i) {
        LOG(DEBUGL) << "Executing path " << i;
        const TGPath &path = paths.getPath(i);

        const Rule &rule = rules[path.ruleid];
        if (rule.getHeads().size() != 1) {
            LOG(ERRORL) << "TGChaseStatic supports only rules with a single"
                " head";
            throw 10;
        }

        //Invoke the execution of the rule using the inputs specified
        TGChase_RuleIO ioRule(path.output, path.inputs);

        std::chrono::system_clock::time_point start =
            std::chrono::system_clock::now();
        //This command will execute a rule, duplicates are not removed
        executeRule(path.ruleid, ioRule);
        std::chrono::duration<double> sec = std::chrono::system_clock::now()
            - start;

        ///Remember the created predicates (later we will need it
        //to clean up the duplicates)
        PredId_t headPredicate = rule.getHeads()[0].getPredicate().getId();
        for(size_t i = 0; i < ioRule.outNodes.size(); ++i) {
            std::string &nameNode = ioRule.outNodes[i];
            size_t sizeNode = ioRule.outSizes[i];
            if (sizeNode > 0) {
                createdPredicates.insert(headPredicate);
                nTriggers += sizeNode;
            }
        }

        StatIteration stat;
        stat.iteration = i;
        stat.rule = &rules[path.ruleid];
        stat.time = sec.count() * 1000;
        stat.derived = false;
        costRules.push_back(stat);
        iterations.insert(std::make_pair(path.output, i));
    }
    LOG(INFOL) << "N. triggers: " << nTriggers;
    stopRun();

#ifdef DEBUG
    //Print some statistics
    std::sort(costRules.begin(), costRules.end());
    int i = 0;
    double sum = 0;
    double sum10 = 0;
    for (auto &el : costRules) {
        LOG(DEBUGL) << "Cost iteration " << el.iteration << " " <<
            el.time;
        i++;
        if (i >= 20)
            break;

        sum += el.time;
        if (i <= 10)
            sum10 += el.time;
    }
    LOG(DEBUGL) << "Sum first 20 rules: " << sum
        << " first 10:" << sum10;
#endif

    LOG(INFOL) << "Removing duplicates ...";
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
    size_t retainedTuples = 0;
    for(auto &predId : createdPredicates) {
        LOG(DEBUGL) << "Removing duplicates for predicate " << predId
            << " " << program->getPredicateName(predId) << " ...";
        retainedTuples += g.mergeNodesWithPredicateIntoOne(predId);
    }
    std::chrono::duration<double> sec = std::chrono::system_clock::now()
        - start;
    LOG(DEBUGL) << "Runtime duplicates removal: "
        << sec.count() * 1000 << "msec";
    LOG(INFOL) << "Retained " << retainedTuples << " tuples";


}
