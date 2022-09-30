#include <vlog/seminaiver.h>
#include <vlog/seminaiver_threaded.h>
#include <vlog/seminaiver_trigger.h>
#include <vlog/wizard.h>
#include <vlog/reasoner.h>
#include <vlog/concepts.h>
#include <vlog/edb.h>
#include <vlog/qsqquery.h>
#include <vlog/qsqr.h>
#include <vlog/exporter.h>

#include <trident/kb/consts.h>
#include <trident/model/table.h>

#include <string>
#include <vector>

long cmpRow(std::vector<uint8_t> *posJoins, const Term_t *row1, const uint64_t *row2) {
    for (int i = 0; i < posJoins->size(); ++i) {
        long r = (row1[i] - row2[(*posJoins)[i]]);
        if (r != 0) {
            return r;
        }
    }
    return 0;
}

void Reasoner::cleanBindings(std::vector<Term_t> &possibleValuesJoins, std::vector<uint8_t> * posJoins, TupleTable *input) {
    if (input != NULL) {
        const size_t sizeBindings = posJoins->size();
        if (possibleValuesJoins.size() / sizeBindings == input->getNRows()) {
            possibleValuesJoins.clear();
        } else {
            std::vector<Term_t> outputBindings;

            //I assume input and bindings are already sorted
            size_t idxInput = 0;
            const uint64_t *currentRowInput = input->getRow(idxInput);
            size_t idxBindings = 0;
            const Term_t *currentRowBindings = &possibleValuesJoins[0];
            while (idxInput < input->getNRows() && idxBindings < possibleValuesJoins.size()) {
                //Compare the two rows
                const long cmpResult = cmpRow(posJoins, currentRowBindings, currentRowInput);
                if (cmpResult < 0) {
                    for (int j = 0; j < sizeBindings; ++j)
                        outputBindings.push_back(currentRowBindings[j]);
                    idxBindings += sizeBindings;
                    if (idxBindings < possibleValuesJoins.size())
                        currentRowBindings = &possibleValuesJoins[idxBindings];
                } else if (cmpResult >= 0) {
                    idxInput++;
                    if (idxInput < input->getNRows())
                        currentRowInput = input->getRow(idxInput);
                }
                if (cmpResult == 0) {
                    idxBindings += sizeBindings;
                    if (idxBindings < possibleValuesJoins.size())
                        currentRowBindings = &possibleValuesJoins[idxBindings];
                }
            }

            //Write all elements that are left
            while (idxBindings < possibleValuesJoins.size()) {
                outputBindings.push_back(possibleValuesJoins[idxBindings++]);
            }
            possibleValuesJoins.swap(outputBindings);
        }
    }
}

size_t Reasoner::estimate(Literal &query, std::vector<uint8_t> *posBindings,
        std::vector<Term_t> *valueBindings, EDBLayer &layer,
        Program &program) {

    QSQQuery rootQuery(query);
    std::unique_ptr<QSQR> evaluator = std::unique_ptr<QSQR>(
            new QSQR(layer, &program));
    TupleTable *cardTable = NULL;
    cardTable = evaluator->evaluateQuery(QSQR_EST, &rootQuery,
            posBindings, valueBindings, true);
    size_t estimate = cardTable->getRow(0)[0];
    delete cardTable;
    return estimate;
}

FCBlock Reasoner::getBlockFromQuery(Literal constantsQuery, Literal &boundQuery,
        std::vector<uint8_t> *posJoins,
        std::vector<Term_t> *possibleValuesJoins) {
    uint8_t nconstants = (uint8_t) constantsQuery.getTupleSize();

    VTuple constantsTuple = constantsQuery.getTuple();

    std::shared_ptr<FCInternalTable> table;
    if (nconstants == 0) {
        table = std::shared_ptr<FCInternalTable>(new SingletonTable(0));
    } else {
        SegmentInserter inserter(nconstants);
        Term_t tuple[256];
        uint8_t nPosToCopy = 0;
        uint8_t posToCopy[256];
        assert(boundQuery.getTupleSize() <= 256);
        for (int i = 0; i < (uint8_t) boundQuery.getTupleSize(); ++i) {
            if (!boundQuery.getTermAtPos(i).isVariable()) {
                posToCopy[nPosToCopy++] = i;
                tuple[i] = boundQuery.getTermAtPos(i).getValue();
            }
        }

        //Add the input tuples
        if (posJoins != NULL) {
            const int addSize = (int) posJoins->size();

            for (size_t i = 0; i < possibleValuesJoins->size(); i += addSize) {
                for (int j = 0; j < addSize; ++j) {
                    tuple[posJoins->at(j)] = possibleValuesJoins->at(i + j);
                }
                inserter.addRow(tuple, posToCopy);
            }
        } else {
            inserter.addRow(tuple, posToCopy);
        }

        std::shared_ptr<const Segment> seg;
        if (inserter.isSorted()) {
            seg = inserter.getSegment();
        } else {
            seg = inserter.getSegment()->sortBy(NULL);
        }
        table =
            std::shared_ptr<FCInternalTable>(
                    new InmemoryFCInternalTable(nconstants, 0, true, seg));

        //change the constantsQuery
        if (possibleValuesJoins != NULL && possibleValuesJoins->size() > 1) {
            uint8_t varId = 1;
            for (int i = 0; i < nconstants; ++i) {
                if (!inserter.getSegment()->getColumn(i)->isConstant()) {
                    constantsTuple.set(VTerm(varId++, 0), i);
                }
            }
        }
    }
    //iteration==1
    return FCBlock(1, table,
            Literal(constantsQuery.getPredicate(), constantsTuple), 0,
            NULL, 0, true);
}

TupleIterator *Reasoner::getIterator(Literal &query,
        std::vector<uint8_t> *posJoins,
        std::vector<Term_t> *possibleValuesJoins,
        EDBLayer &edb, Program &program, bool returnOnlyVars,
        std::vector<uint8_t> *sortByFields) {
    if (posJoins != NULL && possibleValuesJoins != NULL) {
        /* No, let's keep them. --Ceriel
        // Check if there are'nt too many values to check.
        if (possibleValuesJoins->size() > 0) {
        size_t nvals = possibleValuesJoins->size() / posJoins->size();
        if (nvals > 1000) {
        possibleValuesJoins = NULL;
        posJoins = NULL;
        }
        }
        */
    }
    if (query.getPredicate().getType() == EDB) {
        LOG(INFOL) << "Using edb for " << query.tostring(&program, &edb);
        return Reasoner::getEDBIterator(query, posJoins, possibleValuesJoins, edb,
                returnOnlyVars, sortByFields);
    }
    if (posJoins == NULL || posJoins->size() < query.getNVars() || returnOnlyVars || posJoins->size() > 1) {
        ReasoningMode mode = chooseMostEfficientAlgo(query, edb, program, posJoins, possibleValuesJoins);
        if (mode == MAGIC) {
            LOG(INFOL) << "Using magic for " << query.tostring(&program, &edb);
            return Reasoner::getMagicIterator(
                    query, posJoins, possibleValuesJoins, edb, program,
                    returnOnlyVars, sortByFields);
        }
        //top-down
        LOG(INFOL) << "Using top-down for " << query.tostring(&program, &edb);
        return Reasoner::getTopDownIterator(
                query, posJoins, possibleValuesJoins, edb, program,
                returnOnlyVars, sortByFields);
    }

    LOG(INFOL) << "Using incremental reasoning for " << query.tostring(&program, &edb);
    return getIncrReasoningIterator(query, posJoins, possibleValuesJoins, edb, program, returnOnlyVars, sortByFields);
}

TupleIterator *Reasoner::getIncrReasoningIterator(Literal &query,
        std::vector<uint8_t> *posJoins,
        std::vector<Term_t> *possibleValuesJoins,
        EDBLayer &edb, Program &program, bool returnOnlyVars,
        std::vector<uint8_t> *sortByFields) {

    // For now, only works when returnOnlyVars == false;
    assert(! returnOnlyVars);

    // For now, it also only works when posJoins->size() == 1.
    assert(posJoins->size() == 1);
    std::sort(possibleValuesJoins->begin(), possibleValuesJoins->end());

    std::vector<uint8_t> newPosJoins;
    if (posJoins != NULL) {
        newPosJoins = *posJoins;
        for (int j = 0; j < query.getTupleSize(); ++j) {
            if (!query.getTermAtPos(j).isVariable()) {
                //Increment all newPosToJoin
                for (int m = 0; m < newPosJoins.size(); ++m) {
                    if (newPosJoins[m] >= j)
                        newPosJoins[m]++;
                }
            } else {
            }
        }
    }

    //First try an explicit search
    std::string te("TE");
    VTuple t = query.getTuple();
    QSQQuery explQuery(Literal(program.getPredicate(te,
                    Predicate::calculateAdornment(t)), t));


    LOG(DEBUGL) << "Expl query = " << explQuery.getLiteral()->tostring(&program, &edb);

    // edb.query only returns variables, so do the query and add constants.
    TupleTable *tempTable = new TupleTable(query.getNVars());
    edb.query(&explQuery, tempTable, posJoins, possibleValuesJoins);
    TupleTable *tempTable1 = tempTable->sortBy(*posJoins);
    delete tempTable;

    TupleTable *outputTable;
    if (returnOnlyVars) {
        outputTable = tempTable;
    } else {
        outputTable = new TupleTable(3);
        uint64_t val[3];
        val[0] = t.get(0).getValue();
        val[1] = t.get(1).getValue();
        val[2] = t.get(2).getValue();
        for (int idx = 0; idx < tempTable1->getNRows(); idx++) {
            const uint64_t *current = tempTable1->getRow(idx);
            for (int i = 0; i < newPosJoins.size(); i++) {
                val[newPosJoins[i]] = current[i];
            }
        }
        delete tempTable1;
    }


    LOG(DEBUGL) << "Expl query done";

    QSQQuery rootQuery(query);

    //Did I get everything?
    if (outputTable->getNRows() < possibleValuesJoins->size() / posJoins->size()) {
        //Continue with the search
        LOG(DEBUGL) << "Found bindings " << outputTable->getNRows();
        if (outputTable->getNRows() > 0) {
            cleanBindings(*possibleValuesJoins, &newPosJoins, outputTable);
        }

        LOG(INFOL) << "number of possible values = " << (possibleValuesJoins->size() / posJoins->size());
        if (possibleValuesJoins->size() / posJoins->size() <= 1000) {
            // TODO: experiment with threshold on when to use this ...
            //I use QSQR
            QSQR evaluator(edb, &program);
            std::vector<Rule> originalRules = program.getAllRulesByPredicate(
                    rootQuery.getLiteral()->getPredicate().getId());

            //Create a smaller program
            Program clonedProgram = program.clone();
            std::vector<Rule> clonedRules = clonedProgram.getAllRulesByPredicate(
                    rootQuery.getLiteral()->getPredicate().getId());
            //Clean any rule with IDB predicates
            while (clonedRules.size() > 0) {
                if (clonedRules.back().getNIDBPredicates() != 0)
                    clonedRules.pop_back();
                else
                    break;
            }

            LOG(DEBUGL) << "Rules in the cloned program " << clonedRules.size();

            //Execute all simple rules with one IDB in a sequence
            for (std::vector<Rule>::iterator itr = originalRules.begin();
                    itr != originalRules.end() && possibleValuesJoins->size() > 0;
                    ++itr) {

                if (itr->getNIDBPredicates() == 0) {
                    continue;
                } else if (itr->getNIDBPredicates() > 1
                        || possibleValuesJoins->size() == 0) {
                    break;
                }

                //Add the rule
                evaluator.deallocateAllRules();
                LOG(DEBUGL) << "Executing the rule " << itr->tostring(NULL, NULL);
                clonedRules.push_back(*itr);
                evaluator.setProgram(&clonedProgram);

                //Launch only the single rule
                TupleTable *tmpTable = evaluator.evaluateQuery(QSQR_EVAL, &rootQuery, &newPosJoins,
                        possibleValuesJoins, returnOnlyVars);

                if (tmpTable != NULL) {
                    //Clean the bindings
                    TupleTable *tmp1 = tmpTable->sortBy(newPosJoins);
                    delete tmpTable;
                    tmpTable = tmp1;
                    if (tmpTable->getNRows() > 0) {
                        LOG(DEBUGL) << "#values = " << tmpTable->getNRows();
                        cleanBindings(*possibleValuesJoins, &newPosJoins, tmpTable);

                        //Add the temporary bindings to the table
                        tmp1 = outputTable->merge(tmpTable);
                        delete outputTable;
                        delete tmpTable;
                        outputTable = tmp1;
                    }
                }

                //Remove the rule
                clonedRules.pop_back();
            }
        }

        if (possibleValuesJoins->size() > 0) {
            // Decide between magic or QSQ-R, to verify the remaining values.
            int algo = chooseMostEfficientAlgo(query, edb, program,
                    posJoins,
                    possibleValuesJoins);

            TupleIterator *itr = NULL;

            if (algo == TOPDOWN) {
                itr = getTopDownIterator(query, posJoins, possibleValuesJoins,
                        edb, program, returnOnlyVars, &newPosJoins);
            } else {
                itr = getMagicIterator(query, posJoins, possibleValuesJoins,
                        edb, program, returnOnlyVars, &newPosJoins);
            }

            //Add the bindings to a temporary container.
            int rowSize = outputTable->getSizeRow();
            TupleTable *tmp = new TupleTable(rowSize);
            while (itr->hasNext()) {
                itr->next();
                for (int i = 0; i < rowSize; ++i) {
                    tmp->addValue(itr->getElementAt(i));
                }
            }
            //Merge the temporary container into the final one.
            //Keeps the output table sorted.
            TupleTable *tmp1 = outputTable->merge(tmp);
            delete tmp;
            delete outputTable;
            outputTable = tmp1;

            if (itr != NULL)
                delete itr;
        }
    }

    //Return an iterator of the bindings
    std::shared_ptr<TupleTable> pFinalTable(outputTable);

    //Add sort by if requested
    if (sortByFields != NULL && !sortByFields->empty()) {
        std::shared_ptr<TupleTable> sortTab = std::shared_ptr<TupleTable>(
                pFinalTable->sortBy(*sortByFields));
        return new TupleTableItr(sortTab);
    } else {
        return new TupleTableItr(pFinalTable);
    }
}

TupleIterator *Reasoner::getTGMagicIterator(Literal &query,
        EDBLayer &edb, Program &program, bool returnOnlyVars,
        std::string profilerPath,
        std::string storematPath) {

    //To use if the flag returnOnlyVars is set to false
    uint64_t outputTuple[256];    // Used in trident method, so no Term_t
    uint8_t nPosToCopy = 0;
    uint8_t posToCopy[256];
    for (int j = 0; j < query.getTupleSize(); ++j) {
        if (!query.getTermAtPos(j).isVariable()) {
            outputTuple[j] = query.getTermAtPos(j).getValue();
        } else {
            posToCopy[nPosToCopy++] = j;
        }
    }

    //Replace variables with constants if posJoins != NULL.
    VTuple t = query.getTuple();
    VTuple boundTuple = t;
    Predicate pred1(query.getPredicate(), Predicate::calculateAdornment(boundTuple));
    Literal query1(pred1, boundTuple);

    //Get all adorned rules
    std::unique_ptr<Wizard> wizard = std::unique_ptr<Wizard>(new Wizard());
    std::shared_ptr<Program> adornedProgram = wizard->getAdornedProgram(query1, program);
    //Print all rules
#if DEBUG
    LOG(DEBUGL) << "Adorned program:";
    std::vector<Rule> newRules = adornedProgram->getAllRules();
    for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
        LOG(DEBUGL) << itr->tostring(adornedProgram.get(), &edb);
    }
#endif

    //Rewrite and add the rules
    std::pair<PredId_t, PredId_t> inputOutputRelIDs;
    std::shared_ptr<Program> magicProgram = wizard->doMagic(query1, adornedProgram,
            inputOutputRelIDs);

#if DEBUG
    LOG(DEBUGL) << "Magic program:";
    newRules = magicProgram->getAllRules();
    for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
        LOG(DEBUGL) << itr->tostring(magicProgram.get(), &edb);
    }
#endif

    if (profilerPath != "") {
        std::string sout = profilerPath + "-rules";
        std::ofstream out(sout);
        std::vector<Rule> newRules = magicProgram->getAllRules();
        for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
            out << itr->tostring(magicProgram.get(), &edb) << std::endl;
        }
        out.close();
    }

    std::shared_ptr<GBChase> sn = Reasoner::getGBChase(edb, magicProgram.get(), GBChaseAlgorithm::TGCHASE_DYNAMIC_FULLPROV);

    if (profilerPath != "") {
        sn->setPathStoreStatistics(profilerPath);
    }

    //Add all the input tuples in the input relation
    Predicate pred = magicProgram->getPredicate(inputOutputRelIDs.first);
    VTuple onlyConstsTuple(pred.getCardinality());
    int j = 0;
    for (int i = 0; i < t.getSize(); ++i) {
        if (!boundTuple.get(i).isVariable()) {
            onlyConstsTuple.set(VTerm(0, t.get(i).getValue()), j++);
        }
    }

    Literal unboundQuery(pred, onlyConstsTuple);
    sn->getGBGraph().addNode(0, unboundQuery);

    //Exec the materialization
    sn->prepareRun(0, ~0ul);
    sn->run();

    if (storematPath != "") {
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        Exporter exp(sn);
        exp.storeOnFiles(storematPath, true, 0, false);
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        LOG(INFOL) << "Time to index and store the materialization on disk = " << sec.count() << " seconds";
    }

    //Extract the tuples from the output relation
    Literal outputLiteral(magicProgram->getPredicate(inputOutputRelIDs.second), t);
    LOG(DEBUGL) << "outputLiteral = " << outputLiteral.tostring(magicProgram.get(), &edb);
    LOG(DEBUGL) << "returnOnlyVars = " << returnOnlyVars;

    TupleTable *finalTable;
    finalTable = new TupleTable(outputLiteral.getNVars());

    auto itr = sn->getTableItr(outputLiteral.getPredicate().getId());
    std::vector<uint8_t> posVars = outputLiteral.getPosVars();
    while (!itr.isEmpty()) {
        std::shared_ptr<const FCInternalTable> table = itr.getCurrentTable();
        FCInternalTableItr *itrTable = table->getIterator();

        // itrTable contains only variables.
        if (returnOnlyVars) {
            while (itrTable->hasNext()) {
                itrTable->next();
                if (finalTable->getSizeRow() == 0) {
                    Term_t row = 0;
                    finalTable->addRow(&row);
                } else {
                    for (uint8_t j = 0; j < posVars.size(); ++j) {
                        finalTable->addValue(itrTable->getCurrentValue(j));
                    }
                }
            }
        } else {
            while (itrTable->hasNext()) {
                itrTable->next();
                for (uint8_t j = 0; j < nPosToCopy; ++j) {
                    outputTuple[posToCopy[j]] = itrTable->getCurrentValue(j);
                }
                finalTable->addRow(outputTuple);
            }
        }

        table->releaseIterator(itrTable);
        itr.moveNextCount();
    }
    std::shared_ptr<TupleTable> pFinalTable(finalTable);
    return new TupleTableItr(pFinalTable);
}

TupleIterator *Reasoner::getProbMagicIterator(Literal &query,
        EDBLayer &edb, Program &program, bool returnOnlyVars,
        bool optDelProofsStaticAnalysis,
        std::string profilerPath,
        std::string storematPath) {

    //To use if the flag returnOnlyVars is set to false
    uint64_t outputTuple[256];    // Used in trident method, so no Term_t
    uint8_t nPosToCopy = 0;
    uint8_t posToCopy[256];
    for (int j = 0; j < query.getTupleSize(); ++j) {
        if (!query.getTermAtPos(j).isVariable()) {
            outputTuple[j] = query.getTermAtPos(j).getValue();
        } else {
            posToCopy[nPosToCopy++] = j;
        }
    }

    //Replace variables with constants if posJoins != NULL.
    VTuple t = query.getTuple();
    VTuple boundTuple = t;
    Predicate pred1(query.getPredicate(), Predicate::calculateAdornment(boundTuple));
    Literal query1(pred1, boundTuple);

    //Get all adorned rules
    std::unique_ptr<Wizard> wizard = std::unique_ptr<Wizard>(new Wizard());
    std::shared_ptr<Program> adornedProgram = wizard->getAdornedProgram(query1, program);
    //Print all rules
#if DEBUG
    LOG(DEBUGL) << "Adorned program:";
    std::vector<Rule> newRules = adornedProgram->getAllRules();
    for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
        LOG(DEBUGL) << itr->tostring(adornedProgram.get(), &edb);
    }
#endif

    //Rewrite and add the rules
    std::pair<PredId_t, PredId_t> inputOutputRelIDs;
    std::shared_ptr<Program> magicProgram = wizard->doMagic(query1, adornedProgram,
            inputOutputRelIDs);

#if DEBUG
    LOG(DEBUGL) << "Magic program:";
    newRules = magicProgram->getAllRules();
    for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
        LOG(DEBUGL) << itr->tostring(magicProgram.get(), &edb);
    }
#endif

    if (profilerPath != "") {
        std::string sout = profilerPath + "-rules";
        std::ofstream out(sout);
        std::vector<Rule> newRules = magicProgram->getAllRules();
        for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
            out << itr->tostring(magicProgram.get(), &edb) << std::endl;
        }
        out.close();
    }

    std::shared_ptr<GBChase> sn = Reasoner::getProbTGChase(
            edb, magicProgram.get(), optDelProofsStaticAnalysis);

    if (profilerPath != "") {
        sn->setPathStoreStatistics(profilerPath);
    }

    //Add all the input tuples in the input relation
    Predicate pred = magicProgram->getPredicate(inputOutputRelIDs.first);
    VTuple onlyConstsTuple(pred.getCardinality());
    int j = 0;
    for (int i = 0; i < t.getSize(); ++i) {
        if (!boundTuple.get(i).isVariable()) {
            onlyConstsTuple.set(VTerm(0, t.get(i).getValue()), j++);
        }
    }

    Literal unboundQuery(pred, onlyConstsTuple);
    sn->getGBGraph().addNode(0, unboundQuery);

    //Exec the materialization
    sn->prepareRun(0, ~0ul);
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    sn->run();
    std::chrono::duration<double> secMat = std::chrono::system_clock::now() - start;
    LOG(INFOL) << "Runtime materialization = " << secMat.count() * 1000 << " milliseconds";
    LOG(INFOL) << "Derived tuples = " << sn->getNDerivedFacts();
    LOG(INFOL) << "N. nodes = " << sn->getNnodes();
    LOG(INFOL) << "N. edges = " << sn->getNedges();
    LOG(INFOL) << "Triggers = " << sn->getNTriggers();

    if (storematPath != "") {
        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        Exporter exp(sn);
        exp.storeOnFiles(storematPath, true, 0, false);
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        LOG(INFOL) << "Time to index and store the materialization on disk = " << sec.count() << " seconds";
    }

    //Extract the tuples from the output relation
    Literal outputLiteral(magicProgram->getPredicate(inputOutputRelIDs.second), t);
    LOG(DEBUGL) << "outputLiteral = " << outputLiteral.tostring(magicProgram.get(), &edb);
    LOG(DEBUGL) << "returnOnlyVars = " << returnOnlyVars;

    TupleTable *finalTable;
    finalTable = new TupleTable(outputLiteral.getNVars());

    auto itr = sn->getTableItr(outputLiteral.getPredicate().getId());
    std::vector<uint8_t> posVars = outputLiteral.getPosVars();
    while (!itr.isEmpty()) {
        std::shared_ptr<const FCInternalTable> table = itr.getCurrentTable();
        FCInternalTableItr *itrTable = table->getIterator();

        // itrTable contains only variables.
        if (returnOnlyVars) {
            while (itrTable->hasNext()) {
                itrTable->next();
                if (finalTable->getSizeRow() == 0) {
                    Term_t row = 0;
                    finalTable->addRow(&row);
                } else {
                    for (uint8_t j = 0; j < posVars.size(); ++j) {
                        finalTable->addValue(itrTable->getCurrentValue(j));
                    }
                }
            }
        } else {
            while (itrTable->hasNext()) {
                itrTable->next();
                for (uint8_t j = 0; j < nPosToCopy; ++j) {
                    outputTuple[posToCopy[j]] = itrTable->getCurrentValue(j);
                }
                finalTable->addRow(outputTuple);
            }
        }

        table->releaseIterator(itrTable);
        itr.moveNextCount();
    }
    std::shared_ptr<TupleTable> pFinalTable(finalTable);
    return new TupleTableItr(pFinalTable);
}

TupleIterator *Reasoner::getMagicIterator(Literal &query,
        std::vector<uint8_t> *posJoins,
        std::vector<Term_t> *possibleValuesJoins,
        EDBLayer &edb, Program &program, bool returnOnlyVars,
        std::vector<uint8_t> *sortByFields) {


    //To use if the flag returnOnlyVars is set to false
    uint64_t outputTuple[256];    // Used in trident method, so no Term_t
    uint8_t nPosToCopy = 0;
    uint8_t posToCopy[256];
    std::vector<uint8_t> newPosJoins; //This is used because I need the posJoins in the original triple, and not on the variables
    if (posJoins != NULL) {
        newPosJoins = *posJoins;
        for (int j = 0; j < query.getTupleSize(); ++j) {
            if (!query.getTermAtPos(j).isVariable()) {
                //Increment all newPosToJoin
                for (int m = 0; m < newPosJoins.size(); ++m) {
                    if (newPosJoins[m] >= j)
                        newPosJoins[m]++;
                }
                // Line below added (quite important ...) --Ceriel
                outputTuple[j] = query.getTermAtPos(j).getValue();
            } else {
                posToCopy[nPosToCopy++] = j;
            }
        }
    } else {
        for (int j = 0; j < query.getTupleSize(); ++j) {
            if (!query.getTermAtPos(j).isVariable()) {
                outputTuple[j] = query.getTermAtPos(j).getValue();
            } else {
                posToCopy[nPosToCopy++] = j;
            }
        }
    }

    //Replace variables with constants if posJoins != NULL.
    VTuple t = query.getTuple();
    VTuple boundTuple = t;
    if (posJoins != NULL) {
        //The posjoins do not include constants
        int j = 0;
        for (std::vector<uint8_t>::const_iterator itr = newPosJoins.begin();
                itr != newPosJoins.end(); ++itr) {
            boundTuple.set(VTerm(0, possibleValuesJoins->at(j++)), *itr);
        }

    }

    Predicate pred1(query.getPredicate(), Predicate::calculateAdornment(boundTuple));
    Literal query1(pred1, boundTuple);

    //Get all adorned rules
    std::unique_ptr<Wizard> wizard = std::unique_ptr<Wizard>(new Wizard());
    std::shared_ptr<Program> adornedProgram = wizard->getAdornedProgram(query1, program);
    //Print all rules
#if DEBUG
    LOG(DEBUGL) << "Adorned program:";
    std::vector<Rule> newRules = adornedProgram->getAllRules();
    for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
        LOG(DEBUGL) << itr->tostring(adornedProgram.get(), &edb);
    }
#endif

    //Rewrite and add the rules
    std::pair<PredId_t, PredId_t> inputOutputRelIDs;
    std::shared_ptr<Program> magicProgram = wizard->doMagic(query1, adornedProgram,
            inputOutputRelIDs);

#if DEBUG
    LOG(DEBUGL) << "Magic program:";
    newRules = magicProgram->getAllRules();
    for (std::vector<Rule>::iterator itr = newRules.begin(); itr != newRules.end(); ++itr) {
        LOG(DEBUGL) << itr->tostring(magicProgram.get(), &edb);
    }
#endif

    SemiNaiver *naiver = new SemiNaiver(
            edb, magicProgram.get(), true, false, false, -1, false, false) ;

    //Add all the input tuples in the input relation
    Predicate pred = magicProgram->getPredicate(inputOutputRelIDs.first);
    VTuple onlyConstsTuple(pred.getCardinality());
    int j = 0;
    for (int i = 0; i < t.getSize(); ++i) {
        if (!boundTuple.get(i).isVariable()) {
            onlyConstsTuple.set(VTerm(0, t.get(i).getValue()), j++);
        }
    }

    Literal unboundQuery(pred, onlyConstsTuple);
    LOG(DEBUGL) << "unboundQuery = " << unboundQuery.tostring(magicProgram.get(), &edb);
    naiver->addDataToIDBRelation(magicProgram->getPredicate(inputOutputRelIDs.first),
            getBlockFromQuery(unboundQuery, query1,
                newPosJoins.size() != 0 ? &newPosJoins : NULL, possibleValuesJoins));

    //Exec the materialization
    naiver->run(1, 2);

    //Extract the tuples from the output relation
    Literal outputLiteral(magicProgram->getPredicate(inputOutputRelIDs.second), t);
    LOG(DEBUGL) << "outputLiteral = " << outputLiteral.tostring(magicProgram.get(), &edb);
    LOG(DEBUGL) << "returnOnlyVars = " << returnOnlyVars;
    LOG(DEBUGL) << "sortByFields->empty() = " << (sortByFields == NULL ? true : sortByFields->empty());

    FCIterator itr = naiver->getTable(outputLiteral, 0, (size_t) - 1);

    TupleTable *finalTable;
    if (returnOnlyVars) {
        finalTable = new TupleTable(outputLiteral.getNVars());
    } else {
        finalTable = new TupleTable(query.getTupleSize());
    }

    std::vector<uint8_t> posVars = outputLiteral.getPosVars();
    while (!itr.isEmpty()) {
        std::shared_ptr<const FCInternalTable> table = itr.getCurrentTable();
        // LOG(DEBUGL) << "table empty? " << table->isEmpty();
        FCInternalTableItr *itrTable = table->getIterator();

        // itrTable contains only variables.
        if (returnOnlyVars) {
            while (itrTable->hasNext()) {
                itrTable->next();
                if (finalTable->getSizeRow() == 0) {
                    Term_t row = 0;
                    finalTable->addRow(&row);
                } else {
                    for (int j = 0; j < posVars.size(); ++j) {
                        finalTable->addValue(itrTable->getCurrentValue(j));
                    }
                }
                // Not sure about this. Was:
                // for (int j = 0; j < rowsize; ++j) {
                //     finalTable->addValue(itrTable->getCurrentValue(j));
                // }
                // TODO!
            }
        } else {
            while (itrTable->hasNext()) {
                itrTable->next();
                for (int j = 0; j < nPosToCopy; ++j) {
                    outputTuple[posToCopy[j]] = itrTable->getCurrentValue(j);
                }
                finalTable->addRow(outputTuple);
                // LOG(DEBUGL) << "Adding row " << outputTuple[0] << ", " << outputTuple[1] << ", " << outputTuple[2];
            }
        }

        table->releaseIterator(itrTable);
        itr.moveNextCount();
    }

    std::shared_ptr<TupleTable> pFinalTable(finalTable);
    delete naiver;

    if (sortByFields != NULL && !sortByFields->empty()) {
        std::shared_ptr<TupleTable> sortTab = std::shared_ptr<TupleTable>(
                pFinalTable->sortBy(*sortByFields));
        return new TupleTableItr(sortTab);

    } else {
        return new TupleTableItr(pFinalTable);
    }
}

TupleIterator *Reasoner::getMaterializationIterator(Literal &query,
        std::vector<uint8_t> *posJoins,
        std::vector<Term_t> *possibleValuesJoins,
        EDBLayer &edb, Program &program, bool returnOnlyVars,
        std::vector<uint8_t> *sortByFields) {

    Predicate pred = query.getPredicate();
    VTuple tuple = query.getTuple();
    if (pred.getType() == EDB) {
        LOG(INFOL) << "Using edb for " << query.tostring(&program, &edb);
        return Reasoner::getEDBIterator(query, posJoins, possibleValuesJoins, edb,
                returnOnlyVars, sortByFields);
    }

    if (posJoins != NULL) {
        LOG(INFOL) << "getMaterializationIterator with joins not implemented yet";
        throw 10;
    }

    // Run materialization
    SemiNaiver *sn = new SemiNaiver(
            edb, &program, true, false,
            false, -1, false, false);

    sn->run();

    TupleIterator *result = getIteratorWithMaterialization(sn, query, returnOnlyVars, sortByFields);
    delete sn;
    return result;
}

TupleIterator *Reasoner::getIteratorWithMaterialization(SemiNaiver *sn, Literal &query, bool returnOnlyVars,
        std::vector<uint8_t> *sortByFields) {

    FCIterator tableIt = sn->getTableItr(query.getPredicate().getId());
    VTuple tuple = query.getTuple();

    TupleTable *finalTable;
    if (returnOnlyVars) {
        finalTable = new TupleTable(query.getNVars());
    } else {
        finalTable = new TupleTable(query.getTupleSize());
    }

    std::vector<std::pair<uint8_t, uint8_t>> repeated = query.getRepeatedVars();

    while (! tableIt.isEmpty()) {
        std::shared_ptr<const FCInternalTable> table = tableIt.getCurrentTable();
        FCInternalTableItr *itrTable = table->getIterator();
        while (itrTable->hasNext()) {
            itrTable->next();
            bool copy = true;
            for (int i = 0; i < tuple.getSize(); i++) {
                if (! tuple.get(i).isVariable()) {
                    if (itrTable->getCurrentValue(i) != tuple.get(i).getValue()) {
                        copy = false;
                        break;
                    }
                }
            }
            if (copy) {
                for (int i = 0; i < repeated.size(); ++i) {
                    if (itrTable->getCurrentValue(repeated[i].first) != itrTable->getCurrentValue(repeated[i].second)) {
                        copy = false;
                        break;
                    }
                }
            }
            if (! copy) {
                continue;
            }
            if (finalTable->getSizeRow() == 0) {
                Term_t row = 0;
                finalTable->addRow(&row);
            } else {
                for (int i = 0; i < tuple.getSize(); i++) {
                    if (! returnOnlyVars || tuple.get(i).isVariable()) {
                        finalTable->addValue(itrTable->getCurrentValue(i));
                    }
                }
            }
        }
        table->releaseIterator(itrTable);
        tableIt.moveNextCount();
    }

    std::shared_ptr<TupleTable> pFinalTable(finalTable);

    if (sortByFields != NULL && !sortByFields->empty()) {
        std::shared_ptr<TupleTable> sortTab = std::shared_ptr<TupleTable>(
                pFinalTable->sortBy(*sortByFields));
        return new TupleTableItr(sortTab);

    } else {
        return new TupleTableItr(pFinalTable);
    }
}

int Reasoner::getNumberOfIDBPredicates(Literal &query, Program &program) {
    int count = 0;
    std::vector<Rule> rules;
    int idxRules = 0;

    std::unordered_set<Term_t> setQueries;
    std::vector<Literal> queries;
    int idxQueries = 0;
    queries.push_back(query);
    Term_t key = (query.getPredicate().getId() << 16) + query.getPredicate().getAdornment();
    setQueries.insert(key);

    while (idxQueries < queries.size()) {
        Literal lit = queries[idxQueries];

        //Go through all rules and get the ones which match the query
        std::vector<Rule> r = program.getAllRulesByPredicate(lit.getPredicate().getId());
        for (std::vector<Rule>::iterator itr = r.begin(); itr != r.end();
                ++itr) {
            rules.push_back(itr->createAdornment(lit.getPredicate().getAdornment()));
        }

        //Go through all the new rules and get new queries to process
        while (idxRules < rules.size()) {
            Rule *r = &rules[idxRules];
            for (std::vector<Literal>::const_iterator itr = r->getBody().begin();
                    itr != r->getBody().end(); ++itr) {
                Predicate pred = itr->getPredicate();
                if (pred.getType() == IDB) {
                    count++;
                    Term_t key = (pred.getId() << 16) + pred.getAdornment();
                    if (setQueries.find(key) == setQueries.end()) {
                        setQueries.insert(key);
                        queries.push_back(*itr);
                    }
                }
            }
            idxRules++;
        }
        idxQueries++;
    }

    return count;
}

void Reasoner::getMetrics(Literal &query, std::vector<uint8_t> *posBindings, std::vector<Term_t> *valueBindings,
	    EDBLayer &layer, Program &program, Metrics &metrics, int maxDepth, string& idbFeatures) {
    std::unique_ptr<QSQR> evaluator = std::unique_ptr<QSQR>(
            new QSQR(layer, &program));
    memset(&metrics, 0, sizeof(Metrics));
    std::vector<uint32_t> uniqueRules;
    vector<PredId_t> idbIds;
    evaluator->estimateQuery(metrics, maxDepth, query, uniqueRules, idbIds);
    vector<PredId_t> allPredIds =  program.getAllPredicateIDs();
    uint32_t countPreds = allPredIds.size();

    uint8_t card = query.getPredicate().getCardinality();
    bool sub_bound = false;
    bool obj_bound = false;
    for (uint8_t c = 0; c < card; ++c) {
        bool bound = !query.getTermAtPos(c).isVariable();
        if (0 == c) {
            sub_bound = bound;
        } else if (1 == c) {
            obj_bound = bound;
        }
    }
    if (1 == card) {
        obj_bound = sub_bound;
    }

    int i = 0;
    stringstream ss;
    ss << (sub_bound ? "1" : "0");
    ss << ",";
    ss << (obj_bound ? "1" : "0");
    ss << ",";
    for (auto pid : allPredIds) {
        auto it = std::find(idbIds.begin(), idbIds.end(), pid);
        if (it != idbIds.end()) {
            // found
            ss << "1";
        } else {
            ss << "0";
        }
        if (i < countPreds-1){
            ss << ",";
        }
        i += 1;
    }
    idbFeatures = ss.str();
    metrics.countUniqueRules = uniqueRules.size();
}

ReasoningMode Reasoner::chooseMostEfficientAlgo(Literal &query,
        EDBLayer &layer, Program &program,
        std::vector<uint8_t> *posBindings,
        std::vector<Term_t> *valueBindings) {
    uint64_t cost = 0;
    if (posBindings != NULL) {
        //Create a new query with the values substituted
        int idxValues = 0;
        VTuple newTuple = query.getTuple();
        for (std::vector<uint8_t>::iterator itr = posBindings->begin(); itr != posBindings->end();
                ++itr) {
            newTuple.set(VTerm(0, valueBindings->at(idxValues)), *itr);
            idxValues++;
        }
        // Fixed adornments in predicate of literal below.
        Predicate pred1(query.getPredicate(), Predicate::calculateAdornment(newTuple));
        Literal newLiteral(pred1, newTuple);
        size_t singleCost = estimate(newLiteral, NULL, NULL, layer, program);
        LOG(DEBUGL) << "SingleCost is " <<
            singleCost << " nBindings " << (valueBindings->size() / posBindings->size());

        //Are bindings less than 10? Then singleCost is probably about right
        uint64_t nValues = valueBindings->size() / posBindings->size();
        if (nValues > 10) {
            //Copy the first 10 values
            std::vector<Term_t> limitedValueBindings;
            for (int i = 0; i < 10 * posBindings->size(); ++i) {
                limitedValueBindings.push_back(valueBindings->at(i));
            }
            uint64_t tenCost = estimate(query, posBindings,
                    &limitedValueBindings, layer, program);
            LOG(DEBUGL) << "TenCost is " << tenCost;

            //y = mx + b. m is slope, b is constant.
            //I assume the singleCost is the cost at position 0. Otherwise it's not a line
            double m = (double)(tenCost - singleCost) / (10);   // 9? --Ceriel
            long b = singleCost;
            cost = (uint64_t) (m * nValues + b);
        } else {
            cost = singleCost;
        }
    } else {
        cost = estimate(query, NULL, NULL, layer, program);
    }
    ReasoningMode mode = cost < threshold ? TOPDOWN : MAGIC;
    LOG(DEBUGL) << "Deciding whether I should resolve " <<
        query.tostring(&program, &layer) <<
        " with magic or QSQR. Estimated cost: " <<
        cost << " threshold for QSQ-R is " << threshold;
    return mode;
}

TupleIterator *Reasoner::getEDBIterator(Literal &query,
        std::vector<uint8_t> *posJoins,
        std::vector<Term_t> *possibleValuesJoins,
        EDBLayer &edb, bool returnOnlyVars,
        std::vector<uint8_t> *sortByFields) {
    QSQQuery qsqquery(query);
    int nVars = query.getNVars();
    int cardinality = query.getTupleSize();
    TupleTable *table = new TupleTable(nVars);
    edb.query(&qsqquery, table, posJoins, possibleValuesJoins);
    std::shared_ptr<TupleTable> ptable = std::shared_ptr<TupleTable>(table);
    if (! returnOnlyVars && nVars != cardinality) {
        VTuple v = query.getTuple();
        TupleTable *newTable = new TupleTable(cardinality);
        TupleIterator *itr = new TupleTableItr(ptable);
        while (itr->hasNext()) {
            itr->next();
            uint64_t row[256];
            int cnt = 0;
            for (int i = 0; i < cardinality; i++) {
                if (v.get(i).isVariable()) {
                    row[i] = itr->getElementAt(cnt);
                    cnt++;
                } else {
                    row[i] = v.get(i).getValue();
                }
            }
            newTable->addRow(row);
        }
	delete itr;
        ptable = std::shared_ptr<TupleTable>(newTable);
    }
    //Add sort by if requested
    if (sortByFields != NULL && !sortByFields->empty()) {
        std::shared_ptr<TupleTable> sortTab = std::shared_ptr<TupleTable>(
                ptable->sortBy(*sortByFields));
        return new TupleTableItr(sortTab);
    } else {
        return new TupleTableItr(ptable);
    }
}

TupleIterator *Reasoner::getTopDownIterator(Literal &query,
        std::vector<uint8_t> *posJoins,
        std::vector<Term_t> *possibleValuesJoins,
        EDBLayer &edb, Program &program, bool returnOnlyVars,
        std::vector<uint8_t> *sortByFields) {

    LOG(DEBUGL) << "Get topdown iterator for query " << query.tostring(&program, &edb);
    std::vector<uint8_t> newPosJoins;
    if (posJoins != NULL) {
        newPosJoins = *posJoins;
        for (int j = 0; j < query.getTupleSize(); ++j) {
            if (!query.getTermAtPos(j).isVariable()) {
                //Increment all newPosToJoin
                for (int m = 0; m < newPosJoins.size(); ++m) {
                    if (newPosJoins[m] >= j)
                        newPosJoins[m]++;
                }
            } else {
            }
        }
    }

    QSQQuery rootQuery(query);
    LOG(DEBUGL) << "QSQQuery = " << rootQuery.tostring();
    std::unique_ptr<QSQR> evaluator = std::unique_ptr<QSQR>(new QSQR(edb, &program));
    TupleTable *finalTable;
    finalTable = evaluator->evaluateQuery(QSQR_EVAL, &rootQuery, newPosJoins.size() > 0 ? &newPosJoins : NULL,
            possibleValuesJoins, returnOnlyVars);

    //Return an iterator of the bindings
    std::shared_ptr<TupleTable> pFinalTable(finalTable);

    //Add sort by if requested
    if (sortByFields != NULL && !sortByFields->empty()) {
        std::shared_ptr<TupleTable> sortTab = std::shared_ptr<TupleTable>(
                pFinalTable->sortBy(*sortByFields));
        return new TupleTableItr(sortTab);
    } else {
        return new TupleTableItr(pFinalTable);
    }
}

std::shared_ptr<SemiNaiver> Reasoner::getSemiNaiver(EDBLayer &layer,
        Program *p, bool opt_intersect, bool opt_filtering, bool opt_threaded,
        TypeChase typeChase,
        int nthreads, int interRuleThreads, bool shuffleRules, Program *restrictedCheck,
        std::string sameasAlgo) {
    LOG(DEBUGL) << "interRuleThreads = " << interRuleThreads << ", shuffleRules = " << shuffleRules;
    if (interRuleThreads > 0 && restrictedCheck == NULL) {
        std::shared_ptr<SemiNaiver> sn(new SemiNaiverThreaded(
                    layer, p, opt_intersect, opt_filtering,
                    shuffleRules, nthreads, interRuleThreads, sameasAlgo));
        return sn;
    } else {
        std::shared_ptr<SemiNaiver> sn(new SemiNaiver(
                    layer, p, opt_intersect, opt_filtering,
                    opt_threaded, typeChase, nthreads, shuffleRules, false, restrictedCheck, sameasAlgo));
        return sn;
    }
}

std::shared_ptr<TriggerSemiNaiver> Reasoner::getTriggeredSemiNaiver(EDBLayer &layer,
        Program *p,
        TypeChase chase) {
    std::shared_ptr<TriggerSemiNaiver> sn(new TriggerSemiNaiver(
                layer, p, chase));
    return sn;
}

std::shared_ptr<GBChase> Reasoner::getProbTGChase(
        EDBLayer &layer,
        Program *p,
        bool optDelProofsStaticAnalysis) {
    std::shared_ptr<GBChase> sn(new ProbGBChase(layer, *p,
                optDelProofsStaticAnalysis));
    return sn;
}

std::shared_ptr<GBChase> Reasoner::getGBChase(
        EDBLayer &layer,
        Program *p,
        GBChaseAlgorithm typeChase,
        bool queryCont,
        bool edbCheck,
        bool rewriteCliques,
        std::string param1) {

    if ((edbCheck || queryCont) && p->areExistentialRules())
    {
        LOG(WARNL) << "edbCheck and queryCont (query containment) are disabled"
            " because they don't work with existential rules";
        edbCheck = queryCont = false;
    }

    if (typeChase == GBChaseAlgorithm::GBCHASE) {
        std::shared_ptr<GBChase> sn(new GBChase(layer, p, true,
                    GBGraph::ProvenanceType::NOPROV, false));
        return sn;
    } else if (typeChase == GBChaseAlgorithm::TGCHASE_STATIC) {
        std::shared_ptr<GBChase> sn(new TGChaseStatic(layer, p, param1));
        return sn;
    } else if (typeChase == GBChaseAlgorithm::TGCHASE_DYNAMIC) {
        std::shared_ptr<GBChase> sn(new GBChase(layer, p,
                    true, GBGraph::ProvenanceType::NODEPROV,
                    queryCont,
                    edbCheck,
                    rewriteCliques));
        return sn;
    } else if (typeChase == GBChaseAlgorithm::TGCHASE_DYNAMIC_FULLPROV) {
        std::shared_ptr<GBChase> sn(new GBChase(layer, p,
                    true, GBGraph::ProvenanceType::FULLPROV,
                    queryCont,
                    edbCheck,
                    rewriteCliques));
        return sn;
    } else {
        LOG(ERRORL) << "Type of chase is not supported";
        throw 10;
    }
}
