#include <vlog/wizard.h>

#include <kognac/logs.h>

#include <unordered_set>
#include <assert.h>

std::shared_ptr<Program> Wizard::getAdornedProgram(Literal &query, Program &program) {


    std::map<PredId_t, std::vector<uint64_t>> processedQueries;
    std::vector<Literal> queries;
    size_t idxQueries = 0;
    queries.push_back(query);
    std::map<std::pair<uint64_t, uint64_t>, std::vector<Rule>> key2rules;

    while (idxQueries < queries.size()) {
        Literal lit = queries[idxQueries];
        auto predId = lit.getPredicate().getId();
        auto predAdorn = lit.getPredicate().getAdornment();
        bool ok = true;
        if (processedQueries.count(predId)) {
            auto &adornments = processedQueries[predId];
            std::vector<uint64_t> adornmentsToRemove;
            for(auto a : adornments) {
                if (a == predAdorn) {
                    ok = false;
                    break;
                } else {
                    //The following optimation must be double-checked.
                    //I think it's correct, but more testing is needed
                    /*auto card = lit.getPredicate().getCardinality();
                    bool isMoreGeneral = true;
                    bool isMoreSpecific = true;
                    for(size_t i = 0; i < card; ++i) {
                        int qnew = (predAdorn >> i) & 1;
                        int qold = (a >> i) & 1;
                        if (qnew > qold) {
                            //cannot be more general because a has 'b' while
                            //predAdorn has 'f'
                            isMoreGeneral = false;
                        } else if (qnew < qold) {
                            //cannot be more specific because a has 'f' while
                            //predAdorn has 'b'
                            isMoreSpecific = false;
                        } else {
                            //cannot say anything
                        }
                    }
                    assert(isMoreGeneral == false || isMoreSpecific == false);
                    if (isMoreGeneral) {
                        //The new query is more general than an existing one
                        adornmentsToRemove.push_back(a);
                    } else {
                        //There is an adornment that is more general. No need to
                        //process this query
                        ok = false;
                        break;
                    }*/
                }
            }
            for(auto a : adornmentsToRemove) {
                auto k = std::make_pair(predId, a);
                if (key2rules.count(k))
                    key2rules.erase(k);
            }
        }

        if (ok) {
            //Add the query to a cache so that we don't process it anymore
            if (processedQueries.count(predId)) {
                processedQueries[predId].push_back(predAdorn);
            } else {
                std::vector<uint64_t> adornments;
                adornments.push_back(predAdorn);
                processedQueries.insert(std::make_pair(predId, adornments));
            }

            //Go through all rules and get the ones which match the query
            auto rulesIds = program.getRulesIDsByPredicate(lit.getPredicate().getId());
            std::vector<Rule> newRules;
            for (auto ruleId : rulesIds) {
                newRules.push_back(program.getRule(ruleId).
                        createAdornment(lit.getPredicate().getAdornment()));
            }
            key2rules.insert(std::make_pair(std::make_pair(predId, predAdorn), newRules));

            //Go through all the new rules and get new queries to process
            for (Rule &r : newRules) {
                for (std::vector<Literal>::const_iterator itr = r.getBody().begin();
                        itr != r.getBody().end(); ++itr) {
                    Predicate pred = itr->getPredicate();
                    if (pred.getType() == IDB) {
                        queries.push_back(*itr);
                    }
                }
            }
        }
        idxQueries++;
    }

    //Create a new program with the new rules
    std::shared_ptr<Program> newProgram = program.cloneNew();
    newProgram->cleanAllRules();
    std::vector<Rule> rules;
    for(auto &p : key2rules) {
        for(auto rule : p.second) {
            rules.push_back(rule);
        }
    }
    newProgram->addAllRules(rules);
    return newProgram;
}

Literal Wizard::getMagicRelation(const bool hasPriority, std::shared_ptr<Program> newProgram, const Literal &head) {
    const PredId_t pred = head.getPredicate().getId();
    const uint8_t adornment = head.getPredicate().getAdornment();

    //new predicate
    std::string newPred = "MAGICI_" + std::to_string(pred) + "_" + std::to_string(adornment);

    const uint8_t tupleSize = Predicate::getNFields(adornment);
    VTuple newTuple(tupleSize);
    int j = 0;
    uint8_t adornmentNewPred = 0;
    for (uint8_t i = 0; i < head.getTupleSize(); ++i) {
        //Check the "bound" terms
        if ((adornment >> i) & 1) { //Is bound
            adornmentNewPred += 1 << j;
            newTuple.set(head.getTermAtPos(i), j++);
        }
    }
    assert(j == tupleSize);

    const PredId_t newpredid = newProgram->getPredicateID(newPred, tupleSize);
    //LOG(DEBUGL) << "Assigning ID " << (int) newpredid << " to rel " << newPred;
    uint8_t type = IDB;
    if (hasPriority)
        type += 2; //This activates the magic flag in the predicate
    return Literal(Predicate(newpredid, adornmentNewPred, type,
                tupleSize), newTuple);
}

std::shared_ptr<Program> Wizard::doMagic(const Literal &query,
        const std::shared_ptr<Program> inputProgram,
        std::pair<PredId_t, PredId_t> &inputOutputRelIDs) {
    std::shared_ptr<Program> newProgram = inputProgram->cloneNew();
    std::vector<Rule> newRules;
    bool foundPredicate = false;

    std::vector<Rule> originalRules = inputProgram->getAllRules();
    //First pass: add an initial relation to each rule
    for (std::vector<Rule>::iterator itr = originalRules.begin();
            itr != originalRules.end(); ++itr) {

        for (auto head : itr->getHeads()) {
            Literal magicLiteral = getMagicRelation(true, newProgram, head);

            if (head.getPredicate().getId() == query.getPredicate().getId() &&
                    head.getPredicate().getAdornment() == query.getPredicate().getAdornment()) {
                inputOutputRelIDs.first = magicLiteral.getPredicate().getId();
                foundPredicate = true;
            }

            std::vector<Literal> newBody;
            newBody.push_back(magicLiteral);
            //Add all other body literals
            for (std::vector<Literal>::const_iterator itrBody = itr->getBody().begin();
                    itrBody != itr->getBody().end(); ++itrBody) {
                newBody.push_back(*itrBody);
            }
            assert(newBody.size() > 0);
            std::vector<Literal> heads;
            heads.push_back(head);
            Rule r(newRules.size(), heads, newBody, itr->isEGD());
            // LOG(DEBUGL) << "Adding rule " << r.tostring();
            newRules.push_back(r.normalizeVars());
        }
    }

    //Second pass: create an additional rule for each IDB in the rules body
    std::vector<Rule> additionalRules;
    std::vector<std::string> additionalRulesStrings;
    for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end();
            itr++) {
        int npreds = itr->getNIDBPredicates();
        LOG(DEBUGL) << "Processing rule " << itr->tostring();
        assert(npreds > 0);
        if (npreds > 1) { //1 is the one we added in the first step.
            for (int i = 1; i < npreds; ++i) {
                //Put in the body all literals up to the ith IDB literal
                std::vector<Literal> newBody;
                int countIDBs = 0;
                int j = 0;
                for (; j < itr->getBody().size(); ++j) {
                    if (itr->getBody()[j].getPredicate().getType() == IDB) {
                        if (countIDBs == i) {
                            break;
                        }
                        countIDBs++;
                    }
                    newBody.push_back(itr->getBody().at(j));
                }

                //create a new head
                Literal newHead = getMagicRelation(true, newProgram, itr->getBody()[j]);
                assert(newBody.size() > 0);

                if (newBody.size() == 1 &&
                        newBody[0].getPredicate().getId() == newHead.getPredicate().getId() &&
                        newBody[0].getPredicate().getType() == newHead.getPredicate().getType() &&
                        newBody[0].getPredicate().getAdornment() == newHead.getPredicate().getAdornment()) {
                } else {
                    std::vector<Literal> newheads;
                    newheads.push_back(newHead);
                    Rule r(newRules.size() + additionalRules.size(),
                            newheads, newBody, itr->isEGD());
                    Rule normalized_r(r.normalizeVars());
                    std::string s = normalized_r.tostring();
                    bool found = false;
                    for (auto itr : additionalRulesStrings) {
                        if (itr == s) {
                            found = true;
                            break;
                        }
                    }
                    if (! found) {
                        additionalRules.push_back(normalized_r);
                        additionalRulesStrings.push_back(s);
                    } else {
                        LOG(DEBUGL) << "Not adding duplicate rule " << s;
                    }
                    // LOG(DEBUGL) << "Adding rule " << r.tostring();
                }

            }
        }
    }
    std::copy(additionalRules.begin(), additionalRules.end(), std::back_inserter(newRules));

    assert(foundPredicate);
    inputOutputRelIDs.second = query.getPredicate().getId();
    newProgram->cleanAllRules();
    newProgram->addAllRules(newRules);
    return newProgram;
}
