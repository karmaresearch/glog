#include <glog/gbruleexecutor.h>
#include <glog/gbsegment.h>
#include <glog/gbsegmentinserter.h>

std::chrono::duration<double, std::milli> GBRuleExecutor::getDuration(
        DurationType typ) {
    switch (typ) {
        case DUR_FIRST:
            return lastDurationFirst;
        case DUR_MERGE:
            return lastDurationMergeSort;
        case DUR_JOIN:
            return lastDurationJoin;
        case DUR_HEAD:
            return lastDurationCreateHead;
            //case DUR_PREP2TO1:
            //    return lastDurationPrep2to1;
    }
    throw 10;
}

std::string GBRuleExecutor::getStat(StatType typ) {
    switch (typ) {
        case N_BDY_ATOMS:
            return bdyAtoms;
    }
    throw 10;
}

std::shared_ptr<const TGSegment> GBRuleExecutor::projectTuples_structuresharing(
        std::shared_ptr<const TGSegment> tuples,
        const std::vector<int> &posKnownVariables,
        bool isSorted) {
    std::vector<std::shared_ptr<Column>> columns;
    tuples->projectTo(posKnownVariables, columns);
    return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(columns,
                tuples->getNRows(), isSorted, 0, getSegProvenanceType(),
                tuples->getNOffsetColumns()));
}

std::shared_ptr<const TGSegment> GBRuleExecutor::projectTuples(
        std::shared_ptr<const TGSegment> tuples,
        const std::vector<int> &posKnownVariables) {
    auto pc = tuples->getProvenanceType();
    assert(provenanceType != GBGraph::ProvenanceType::FULLPROV || pc == SEG_FULLPROV);

    if (posKnownVariables.size() == 1) {
        if (shouldTrackProvenance()) {
            assert(tuples->getProvenanceType() != SEG_NOPROV);
            if (tuples->getProvenanceType() == SEG_SAMENODE) {
                std::vector<Term_t> projection;
                tuples->projectTo(posKnownVariables[0], projection);
                return std::shared_ptr<const TGSegment>(
                        new UnaryWithConstProvTGSegment(
                            projection, tuples->getNodeId(), false, 0));

            } else {
                std::vector<std::pair<Term_t, Term_t>> projection;
                tuples->projectTo(posKnownVariables[0],
                        projection);
                return std::shared_ptr<const TGSegment>(
                        new UnaryWithProvTGSegment(
                            projection, ~0ul, false, 0));
            }
        } else {
            std::vector<Term_t> projection;
            tuples->projectTo(posKnownVariables[0], projection);
            return std::shared_ptr<const TGSegment>(new UnaryTGSegment(
                        projection, ~0ul, false, 0));
        }
    } else if (posKnownVariables.size() == 2) {
        if (shouldTrackProvenance()) {
            assert(tuples->getProvenanceType() != SEG_NOPROV);
            if (tuples->getProvenanceType() == SEG_SAMENODE) {
                std::vector<std::pair<Term_t, Term_t>> projection;
                tuples->projectTo(posKnownVariables[0],
                        posKnownVariables[1],
                        projection);
                return std::shared_ptr<const TGSegment>(
                        new BinaryWithConstProvTGSegment(
                            projection, tuples->getNodeId(), false, 0));
            } else {
                std::vector<BinWithProv> projection;
                tuples->projectTo(posKnownVariables[0],
                        posKnownVariables[1],
                        projection);
                return std::shared_ptr<const TGSegment>(
                        new BinaryWithProvTGSegment(
                            projection, ~0ul, false, 0));
            }
        } else {
            std::vector<std::pair<Term_t, Term_t>> projection;
            tuples->projectTo(posKnownVariables[0],
                    posKnownVariables[1],
                    projection);
            return std::shared_ptr<const TGSegment>(new BinaryTGSegment(
                        projection, ~0ul, false, 0));
        }
    } else { //0 or more than 2
        std::vector<std::shared_ptr<Column>> columns;
        tuples->projectTo(posKnownVariables, columns);
        bool remainSorted = tuples->isSorted();
        for(int i = 0; i < posKnownVariables.size(); ++i) {
            if (posKnownVariables[i] != i) {
                remainSorted = false;
                break;
            }
        }
        //auto npc = getSegProvenanceType();
        //assert(npc == pc);
        return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(columns,
                    tuples->getNRows(), remainSorted, 0,
                    pc, tuples->getNOffsetColumns()));
    }
}

std::shared_ptr<const TGSegment> GBRuleExecutor::projectHead(
        const Literal &head,
        std::vector<size_t> &vars,
        std::shared_ptr<const TGSegment> tuples,
        bool shouldSort,
        bool shouldDelDupl,
        std::vector<std::shared_ptr<Column>> &intermediateResultsNodes,
        GBRuleInput &node) {
    const auto &tupleHead = head.getTuple();
    std::vector<int> posKnownVariables;
    std::vector<std::pair<size_t, Term_t>> posConstants;
    for(int i = 0; i < tupleHead.getSize(); ++i) {
        auto varId = tupleHead.get(i).getId();
        bool found = false;
        for(int j = 0; j < vars.size(); ++j)
            if (vars[j] == varId) {
                found = true;
                posKnownVariables.push_back(j);
                break;
            }
        if (!found) {
            if (!tupleHead.get(i).isVariable()) {
                posConstants.push_back(std::make_pair(
                            i, tupleHead.get(i).getValue()));
            } else {
                LOG(ERRORL) << "Should not happen. At this stage, all variables "
                    "are known";
                throw 10;
            }
        }
    }

    if (posKnownVariables.size() == 0) {
        tuples = projectTuples(tuples, posKnownVariables);
    } else if (posKnownVariables.size() == 1) {
        if (tuples->getNColumns() == 1) {
            //Do nothing
        } else {
            //Must do a unary projection
            tuples = projectTuples(tuples, posKnownVariables);
        }
    } else if (posKnownVariables.size() == 2) {
        if (tuples->getNColumns() == 2) {
            //I might have to swap the columns
            if (posKnownVariables[0] == 1 && posKnownVariables[1] == 0) {
                tuples = tuples->swap();
            } else {
                //Do nothing
                assert(posKnownVariables[0] == 0 && posKnownVariables[1] == 1);
            }
        } else {
            //Must do a binary projection
            tuples = projectTuples(tuples, posKnownVariables);
        }
    } else {
        tuples = projectTuples(tuples, posKnownVariables);
    }

    //Clean up the duplicates if necessary
    if (shouldSort) {
        if (shouldDelDupl) {
            tuples = tuples->sort();
            tuples = tuples->unique();
        } else {
            tuples = tuples->sort();
        }
    } else {
        if (shouldDelDupl) {
            assert(intermediateResultsNodes.size() == 0);
            tuples = tuples->unique();
        }
    }

    if (!posConstants.empty()) {
        tuples = addConstants(tuples, posConstants);
    }

    return tuples;
}

std::shared_ptr<const TGSegment> GBRuleExecutor::addConstants(
        std::shared_ptr<const TGSegment> tg,
        std::vector<std::pair<size_t, Term_t>> constants)
{
    size_t nTerms = tg->getNColumns() + constants.size();
    size_t sizeRow = nTerms + tg->getNOffsetColumns();
    std::unique_ptr<Term_t[]> row = std::unique_ptr<Term_t[]>(
            new Term_t[sizeRow]);
    std::vector<std::pair<size_t, size_t>> posVars;

    size_t literalSize = tg->getNColumns() + constants.size();
    size_t idxPosConsts = 0;
    size_t idxPosVars = 0;
    for(size_t i = 0; i < literalSize; ++i) {
        if (idxPosConsts < constants.size() &&
                i == constants[idxPosConsts].first) {
            row[i] = constants[idxPosConsts].second;
            idxPosConsts++;
        } else if (idxPosVars < tg->getNColumns()) {
            posVars.push_back(std::make_pair(idxPosVars, i));
            idxPosVars++;
        }
    }
    assert(constants.size() == idxPosConsts);
    assert(idxPosVars == tg->getNColumns());
    std::unique_ptr<GBSegmentInserter> ins = GBSegmentInserter::getInserter(
            sizeRow, tg->getNOffsetColumns(), false);
    auto itr = tg->iterator();
    while (itr->hasNext()) {
        itr->next();
        //Copy the variables
        for(auto &p : posVars) {
            row[p.second] = itr->get(p.first);
        }
        if (tg->getNOffsetColumns() > 0) {
            row[nTerms] = itr->getNodeId();
            for(size_t i = 1; i < tg->getNOffsetColumns(); ++i) {
                row[nTerms + i] = itr->getProvenanceOffset(i-1);
            }
        }
        ins->add(row.get());
    }

    return ins->getSegment(tg->getNodeId(),
            tg->isSorted(), 0, tg->getProvenanceType(),
            tg->getNOffsetColumns());
}


void GBRuleExecutor::computeVarPos(std::vector<size_t> &leftVars,
        int bodyAtomIdx,
        const std::vector<Literal> &bodyAtoms,
        const std::vector<Literal> &heads,
        std::vector<std::pair<int, int>> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight) {
    std::set<int> futureOccurrences;
    for(size_t i = bodyAtomIdx + 1; i < bodyAtoms.size(); ++i) {
        auto &lit = bodyAtoms[i];
        auto vars = lit.getAllVars();
        for(auto v : vars)
            futureOccurrences.insert(v);
    }
    for(const auto &head : heads) {
        auto headVars = head.getAllVars();
        for(auto v : headVars)
            futureOccurrences.insert(v);
    }

    auto &rightBodyAtom = bodyAtoms[bodyAtomIdx];
    auto rightVars = rightBodyAtom.getAllVarsAndPos();
    std::set<int> addedVars;
    if (bodyAtomIdx > 0) {
        //First check whether variables in the left are used in the next body
        //atoms or in the head
        for(size_t i = 0; i < leftVars.size(); ++i) {
            if (futureOccurrences.count(leftVars[i])) {
                copyVarPosLeft.push_back(i);
                addedVars.insert(leftVars[i]);
            }
        }
        //Search for the join variable
        bool found = false;
        for(int i = 0; i < rightVars.size(); ++i) {
            for(int j = 0; j < leftVars.size(); ++j) {
                if (rightVars[i].first == leftVars[j]) {
                    joinVarPos.push_back(std::pair<int,int>(j,
                                rightVars[i].second));
                    found = true;
                    break;
                }
            }
        }
    }

    //Copy variables from the right atom
    for(size_t i = 0; i < rightVars.size(); ++i) {
        if (futureOccurrences.count(rightVars[i].first) &&
                !addedVars.count(rightVars[i].first)) {
            copyVarPosRight.push_back(rightVars[i].second);
            addedVars.insert(rightVars[i].first);
        }
    }
}

void GBRuleExecutor::addBuiltinFunctions(
        std::vector<BuiltinFunction> &out,
        const std::vector<Literal> &atoms,
        const std::vector<size_t> &vars) {
    for(auto &atom : atoms) {
        auto fn = layer.getBuiltinFunction(atom);
        bool allVarsAreFound = true;
        int idx = 0;
        for(int i = 0; i < atom.getTupleSize(); ++i) {
            auto term = atom.getTermAtPos(i);
            if (term.isVariable()) {
                //Search among vars
                bool found = false;
                for(int j = 0; j < vars.size(); ++j) {
                    if (term.getId() == vars[j]) {
                        fn.posArgs[idx++] = j;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    allVarsAreFound = false;
                    break;
                }
            }
        }
        if (allVarsAreFound) {
            out.push_back(fn);
        }
    }
}

void GBRuleExecutor::shouldSortAndRetainEDBSegments(
        bool &shouldSort,
        bool &shouldRetainUnique,
        const Literal &atom,
        std::vector<int> &copyVarPos) {
    if (copyVarPos.size() == 0) {
        shouldSort = false;
        shouldRetainUnique = false;
        return;
    }

    //Sort should be done if the order of variables does not reflect the one of
    //the literal (we assume that the relation is sorted left-to-right

    //Watch out, this procedure does not work with repeated variables
    int nconsts = 0;
    int nvars = 0;
    int idxCopyVar = 0;
    shouldSort = false;
    for(size_t i = 0; i < atom.getTupleSize(); ++i) {
        auto t = atom.getTermAtPos(i);
        if (t.isVariable()) {
            nvars++;
            if (idxCopyVar < copyVarPos.size()) {
                if (copyVarPos[idxCopyVar] == i) {
                    idxCopyVar++;
                } else {
                    shouldSort = true;
                }
            }
        } else {
            nconsts++;
        }
    }
    shouldRetainUnique = copyVarPos.size() < nvars;
}

std::shared_ptr<const TGSegment> GBRuleExecutor::processAtom_IDB(
        const Literal &atom, std::vector<size_t> &nodeIdxs,
        std::vector<int> &copyVarPos,
        bool lazyMode,
        bool replaceOffsets) {
    if (atom.getNConstants() > 0) {
        std::vector<Term_t> filterConstants;
        for(size_t i = 0; i < atom.getTupleSize(); ++i) {
            if (atom.getTermAtPos(i).isVariable()) {
                filterConstants.push_back(~0ul);
            } else {
                filterConstants.push_back(atom.getTermAtPos(i).getValue());
            }
        }
        auto intermediateResults = g.mergeNodes(nodeIdxs, filterConstants,
                copyVarPos, lazyMode, replaceOffsets);
        return intermediateResults;
    } else {
        auto intermediateResults = g.mergeNodes(
                nodeIdxs, copyVarPos, lazyMode, replaceOffsets);
        return intermediateResults;
    }
}

std::shared_ptr<const TGSegment> GBRuleExecutor::processAtom_EDB(
        const Literal &atom,
        std::vector<int> &copyVarPos) {
    PredId_t p = atom.getPredicate().getId();
    if (!edbTables.count(p)) {
        auto edbTable = layer.getEDBTable(p);
        edbTables[p] = edbTable;
    }

    if (atom.getTupleSize() == 0) {
        LOG(ERRORL) << "Atoms with arity 0 are not supported";
        throw 10;
    }

    if (atom.hasRepeatedVars()) {
        LOG(ERRORL) << "For now we do not handle atoms with repeated vars";
        throw 10;
    }

    auto &table = edbTables[p];

    //I can apply structure sharing safetly
    bool shouldSort = true;
    bool shouldRetainUnique = true;
    std::vector<std::shared_ptr<Column>> columns;
    size_t nrows = ~0ul;
    if (table->useSegments()) {
        if (atom.getNConstants() == 0)
        {
            auto seg = table->getSegment();
            nrows = seg->getNRows();
            for(int i = 0; i < atom.getTupleSize(); ++i) {
                auto col = seg->getColumn(i);
                columns.push_back(col);
            }
        } else {
            std::vector<Term_t> constants;
            std::vector<size_t> posConstants;
            for(size_t i = 0; i < atom.getTupleSize(); ++i)
            {
                auto t = atom.getTermAtPos(i);
                if (!t.isVariable())
                {
                    constants.push_back(t.getValue());
                    posConstants.push_back(i);
                }
            }
            auto seg = table->getSegment()->filter(constants, posConstants);
            nrows = seg->getNRows();
            for(int i = 0; i < atom.getTupleSize(); ++i) {
                auto col = seg->getColumn(i);
                columns.push_back(col);
            }
        }
    } else {
        //e.g., Trident database
        nrows = table->getCardinality(atom);
        std::vector<uint8_t> presortPos;
        int varIdx = 0;
        for(size_t i = 0; i < atom.getTupleSize(); ++i) {
            auto term = atom.getTermAtPos(i);
            std::shared_ptr<Column> col;
            if (term.isVariable()) {
                col = std::shared_ptr<Column>(
                        new EDBColumn(layer, atom, varIdx, presortPos,
                            false));
                presortPos.push_back(varIdx);
                varIdx++;
            } else {
                col = std::shared_ptr<Column>(
                        new CompressedColumn(term.getValue(), nrows));
            }
            columns.push_back(col);
        }
    }

    shouldSortAndRetainEDBSegments(
            shouldSort,
            shouldRetainUnique,
            atom,
            copyVarPos);

    size_t nProvenanceColumns = 0;
    if (shouldTrackProvenance()) {
        columns.push_back(std::shared_ptr<Column>(
                    new CompressedColumn(~0ul, nrows)));
        nProvenanceColumns = 1;
        if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
            columns.push_back(std::shared_ptr<Column>(
                        new CompressedColumn(0, nrows, 1)));
            nProvenanceColumns = 2;
        }
    }

    auto seg = std::shared_ptr<const TGSegment>(
            new TGSegmentLegacy(columns, nrows, true, 0,
                getSegProvenanceType(), nProvenanceColumns));
    std::shared_ptr<const TGSegment> output;

    if (copyVarPos.size() != atom.getTupleSize()) {
        //There is a projection. We might have to remove duplicates
        auto projectedSegment = projectTuples_structuresharing(seg, copyVarPos,
                !shouldSort);
        std::shared_ptr<const TGSegment> sortedSegment;
        if (shouldSort) {
            sortedSegment = projectedSegment->sort();
        } else {
            sortedSegment = projectedSegment;
        }
        std::shared_ptr<const TGSegment> uniqueSegment;
        if (shouldRetainUnique) {
            uniqueSegment = sortedSegment->unique();
        } else {
            uniqueSegment = sortedSegment;
        }
        output = uniqueSegment;
    } else {
        output = seg;
    }

    if (loadAllEDB) {
        size_t ncols = output->getNColumns() + output->getNOffsetColumns();
        auto inserter = GBSegmentInserter::getInserter(ncols,
                output->getNOffsetColumns(), false);
        std::unique_ptr<Term_t[]> row = std::unique_ptr<Term_t[]>(new Term_t[ncols]);
        auto itr = output->iterator();
        while (itr->hasNext()) {
            itr->next();
            for(size_t i = 0; i < output->getNColumns(); ++i) {
                row[i] = itr->get(i);
            }
            row[output->getNColumns()] = itr->getNodeId();
            for(size_t i = 1; i < output->getNOffsetColumns(); ++i) {
                row[output->getNColumns() + i] = itr->getProvenanceOffset(i-1);
            }
            inserter->add(row.get());
        }
        output = inserter->getSegment(output->getNodeId(), output->isSorted(), 0,
                output->getProvenanceType(), output->getNOffsetColumns());
    }

    return output;
}

void GBRuleExecutor::shouldSortDelDupls(const Literal &head,
        const std::vector<Literal> &bodyAtoms,
        const std::vector<std::vector<size_t>> &bodyNodes,
        bool &shouldSort,
        bool &shouldDelDupl) {
    if (bodyAtoms.size() == 1) {
        auto &bodyAtom = bodyAtoms[0];
        auto th = head.getTuple();
        auto tb = bodyAtom.getTuple();

        std::set<uint32_t> varsInHead;
        for(int i = 0; i < th.getSize(); ++i) {
            if (th.get(i).isVariable())
                varsInHead.insert(th.get(i).getId());
        }
        std::vector<uint32_t> listVarsInBody;
        std::set<uint32_t> varsInBody;
        for(int i = 0; i < tb.getSize(); ++i) {
            if (tb.get(i).isVariable()) {
                if (!varsInBody.count(tb.get(i).getId())) {
                    varsInBody.insert(tb.get(i).getId());
                    listVarsInBody.push_back(tb.get(i).getId());
                }
            }
        }
        size_t nsharedvars = 0;
        for(auto v : varsInBody)
            if (varsInHead.count(v))
                nsharedvars++;
        shouldDelDupl = nsharedvars != varsInBody.size();

        uint32_t prevVar = -1;
        size_t posInBody = 0;
        int i = 0;
        bool sortOk = true;
        for(; i < th.getSize(); ++i) {
            if (th.get(i).isVariable() && th.get(i).getId() != prevVar) {
                if (varsInBody.count(th.get(i).getId())) {
                    if (listVarsInBody[posInBody] !=
                            th.get(i).getId()) {
                        sortOk = false;
                        break;
                    } else {
                        posInBody++;
                        if (posInBody == listVarsInBody.size())
                            break;
                    }
                    prevVar = th.get(i).getId();
                } else {
                    //There is an existential variable. Force sorting...
                    sortOk = false;
                    break;
                }
            }
        }
        //This check verifies whether there are some trailing existential
        //variables (which require sorting)
        sortOk = sortOk && i == th.getSize();
        shouldSort = (!sortOk) ||
            (bodyNodes.size() > 0 && bodyNodes[0].size() > 1);
        if (!retainUnique)
            shouldDelDupl = false;
    } else {
        shouldSort = true;
        if (retainUnique)
            shouldDelDupl = true;
        else
            shouldDelDupl = false;
    }
}

bool __sorter_extvar_dep(
        const std::pair<uint8_t, std::vector<uint8_t>> &a,
        const std::pair<uint8_t, std::vector<uint8_t>> &b) {
    auto minLen = std::min(a.second.size(), b.second.size());
    for(int i = 0; i < minLen; ++i) {
        if (a.second[i] < b.second[i])
            return true;
        else if (a.second[i] > b.second[i]) {
            return false;
        }
    }
    return a.second.size() > b.second.size();
}

bool isprefix(const std::vector<uint8_t> &a,
        const std::vector<uint8_t> &b) {
    if (a.size() > b.size())
        return false;
    for(int i = 0; i < a.size(); ++i) {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

std::shared_ptr<const TGSegment> GBRuleExecutor::addExistentialVariables(
        Rule &rule,
        std::shared_ptr<const TGSegment> tuples,
        std::vector<size_t> &vars) {
    assert(tuples->getNRows() > 0);
    assert(provenanceType != GBGraph::ProvenanceType::FULLPROV); // Not supported

    //Get list existential variables
    std::set<size_t> extvars;
    for(auto &v :rule.getExistentialVariables())
        extvars.insert(v);

    //Assign null values to each variable
    assert(extvars.size() > 0);
    //First compute the frontier variables for each existential variable
    std::map<uint8_t, std::vector<uint8_t>> depVars;
    for(auto &headAtom : rule.getHeads()) {
        std::set<uint8_t> frontierAtomVars;
        std::set<uint8_t> extAtomVars;
        for (auto var : headAtom.getAllVars()) {
            if (extvars.count(var)) {
                extAtomVars.insert(var);
            } else {
                frontierAtomVars.insert(var);
            }
            for(auto var : extAtomVars) {
                if (!depVars.count(var)) {
                    depVars.insert(std::make_pair(var, std::vector<uint8_t>()));
                }
                std::copy(frontierAtomVars.begin(), frontierAtomVars.end(),
                        std::back_inserter(depVars[var]));
            }
        }
    }
    std::vector<std::pair<uint8_t, std::vector<uint8_t>>> newDepVars;
    for (auto &entry : depVars) {
        std::sort(entry.second.begin(), entry.second.end());
        auto last = std::unique(entry.second.begin(), entry.second.end());
        entry.second.erase(last, entry.second.end());
        newDepVars.push_back(std::make_pair(entry.first, entry.second));
    }
    //Sort the dependencies to optimize the sorting
    std::sort(newDepVars.begin(), newDepVars.end(), __sorter_extvar_dep);
    const auto nrows = tuples->getNRows();

    //For each existential variable, assign null values.
    //The other values are copied
    auto counter = g.getCounterNullValues() + 1;
    for(int i = 0; i < newDepVars.size(); ++i) {
        vars.push_back(newDepVars[i].first);
        std::vector<uint8_t> depVarPos;
        for(int j = 0; j < newDepVars[i].second.size(); ++j) {
            for(int m = 0; m < vars.size(); ++m) {
                if (vars[m] == newDepVars[i].second[j]) {
                    depVarPos.push_back(m);
                    break;
                }
            }
        }
        if (i == 0 || !isprefix(newDepVars[i].second,
                    newDepVars[i - 1].second)) {
            tuples = tuples->sortBy(depVarPos);
        }

        //These are three possible containers for the final result. They will
        //contain the list of known variables and the existentially quantified
        //ones
        std::vector<Term_t> terms1;
        std::vector<std::pair<Term_t, Term_t>> terms2;
        std::vector<BinWithProv> terms3;
        std::unique_ptr<SegmentInserter> terms4;
        std::unique_ptr<Term_t> terms4_row;
        const int nfields = tuples->getNColumns() + 1;
        int mode;
        if (nfields == 1) {
            mode = !shouldTrackProvenance() ? 0 : tuples->getProvenanceType() == SEG_SAMENODE ?
                1 : 2;
        } else if (nfields == 2) {
            mode = !shouldTrackProvenance() ? 3 : tuples->getProvenanceType() == SEG_SAMENODE ?
                4 : 5;
        } else {
            if (!shouldTrackProvenance()) {
                mode = 6;
                terms4 = std::unique_ptr<SegmentInserter>(new
                        SegmentInserter(nfields));
                terms4_row = std::unique_ptr<Term_t>(new Term_t[nfields]);
            } else {
                mode = 7;
                terms4 = std::unique_ptr<SegmentInserter>(new
                        SegmentInserter(nfields + 1));
                terms4_row = std::unique_ptr<Term_t>(new Term_t[nfields + 1]);
            }
        }

        switch(mode) {
            case 1:
            case 2:
                terms1.resize(nrows);
                break;
            case 3:
            case 4:
                terms2.resize(nrows);
                break;
            case 5:
                terms3.resize(nrows);
                break;
        }

        const auto nvars = depVarPos.size();
        std::vector<Term_t> prevTerms(nvars, ~0ul);
        auto itr = tuples->iterator();
        size_t currentRow = 0;
        while (itr->hasNext()) {
            itr->next();
            //Compute a suitable null value for the existentially quantified
            //variable
            bool equalPrev = true;
            for(int j = 0; j < depVarPos.size(); ++j) {
                const auto varPos = depVarPos[j];
                if (itr->get(varPos) != prevTerms[j]) {
                    equalPrev = false;
                    break;
                }
            }
            if (!equalPrev) {
                counter++;
                for(int j = 0; j < depVarPos.size(); ++j) {
                    const auto varPos = depVarPos[j];
                    prevTerms[j] = itr->get(varPos);
                }
            }
            //Copy the previous values (and provenance) and the existential
            //value
            switch (mode) {
                case 0:
                case 1:
                    terms1[currentRow] = counter; break;
                case 2:
                    terms2[currentRow].second = itr->getNodeId();
                case 3:
                case 4:
                    if (tuples->getNColumns() == 0) {
                        terms2[currentRow].first = counter;
                    } else {
                        terms2[currentRow].first = itr->get(0);
                        terms2[currentRow].second = counter;
                    }
                    break;
                case 5:
                    terms3[currentRow].first = itr->get(0);
                    terms3[currentRow].second = counter;
                    terms3[currentRow].node = itr->getNodeId();
                    break;
                case 6:
                case 7:
                    for(int i = 0; i < tuples->getNColumns(); ++i) {
                        terms4_row.get()[i] = itr->get(i);
                    }
                    terms4_row.get()[nfields - 1] = counter;
                    if (mode == 7) {
                        terms4_row.get()[nfields] = itr->getNodeId();
                    }
                    terms4->addRow(terms4_row.get());
                    break;
            }
            currentRow++;
        }
        //Create new tuples object
        switch (mode) {
            case 0:
                tuples = std::shared_ptr<const TGSegment>(
                        new UnaryTGSegment(terms1, ~0ul, false, 0));
                break;
            case 1:
                tuples = std::shared_ptr<const TGSegment>(
                        new UnaryWithConstProvTGSegment(terms1,
                            tuples->getNodeId(), false, 0));
                break;
            case 2:
                tuples = std::shared_ptr<const TGSegment>(
                        new UnaryWithProvTGSegment(terms2, ~0ul, false, 0));
                break;

            case 3:
                tuples = std::shared_ptr<const TGSegment>(
                        new BinaryTGSegment(terms2, ~0ul, false, 0));
                break;
            case 4:
                tuples = std::shared_ptr<const TGSegment>(
                        new BinaryWithConstProvTGSegment(terms2,
                            tuples->getNodeId(), false, 0));
                break;
            case 5:
                tuples = std::shared_ptr<const TGSegment>(
                        new BinaryWithProvTGSegment(terms3, ~0ul, false, 0));
                break;
            case 6:
            case 7:
                auto seg = terms4->getSegment();
                std::vector<std::shared_ptr<Column>> columns;
                size_t offsetColumns = (shouldTrackProvenance()) ? 1 : 0;
                auto nfieldsToCopy = nfields + offsetColumns;
                for(int i = 0; i < nfieldsToCopy; ++i) {
                    columns.push_back(seg->getColumn(i));
                }
                tuples = std::shared_ptr<const TGSegment>(
                        new TGSegmentLegacy(columns, seg->getNRows(), false,
                            0, getSegProvenanceType(), offsetColumns));
                break;
        }
    }
    g.setCounterNullValues(counter);
    return tuples;
}

std::shared_ptr<const TGSegment> GBRuleExecutor::performRestrictedCheck(
        Rule &rule,
        const std::shared_ptr<const TGSegment> tuples,
        const std::vector<size_t> &varTuples) {
    assert(provenanceType != GBGraph::ProvenanceType::FULLPROV); //not supported yet
    std::shared_ptr<const TGSegment> newtuples = tuples;
    for(auto &headAtom : rule.getHeads()) {
        if (tuples->isEmpty())
            break;

        std::vector<size_t> nodes; //Not needed (no caching)
        //Find join variables
        std::vector<std::pair<int, int>> joinVarPos;
        std::vector<int> varsToCopyRight;
        for(int  i = 0; i < headAtom.getTupleSize(); ++i) {
            auto t = headAtom.getTermAtPos(i);
            if (t.isVariable()) {
                for(int j = 0; j < varTuples.size(); ++j) {
                    if (t.getId() == varTuples[j]) {
                        varsToCopyRight.push_back(i);
                        joinVarPos.push_back(std::make_pair(
                                    j,joinVarPos.size()));
                        break;
                    }
                }
            } else {
                LOG(ERRORL) << "Not sure what will happen if there are"
                    " constants in the head during the restriction check"
                    ". Throw an exception ...";
                throw 10;
            }
        }

        //We want to retain all variables
        std::vector<int> copyVarPosLeft;
        for(int j = 0; j < varTuples.size(); ++j) {
            copyVarPosLeft.push_back(j);
        }

        if (g.areNodesWithPredicate(headAtom.getPredicate().getId())) {
            const auto &nodesRight = g.getNodeIDsWithPredicate(
                    headAtom.getPredicate().getId());
            std::shared_ptr<const TGSegment> inputRight = g.mergeNodes(
                    nodesRight,
                    varsToCopyRight);

            //Prepare the container that will store the retained tuples
            const int extraColumns = shouldTrackProvenance() ? 1 : 0;
            int nfields = tuples->getNColumns() + extraColumns;
            std::unique_ptr<GBSegmentInserter> outputJoin = GBSegmentInserter::
                getInserter(nfields, extraColumns, false);

            if (tuples->getNColumns() == 0) {
                //Frontier is empty. In this case, if there are non-empty
                //nodes (as it is), then the check fails immediately.
                assert(copyVarPosLeft.size() == 0);
                //Create an empty set of tuples
            } else {
                //Perform a left join
                leftjoin(false,
                        false,
                        tuples,
                        nodes,
                        inputRight,
                        joinVarPos,
                        copyVarPosLeft,
                        outputJoin,
                        true);
            }

            //Create a TGSegment from SegmentInserter
            //std::shared_ptr<const Segment> seg = outputJoin->getSegment();
            //tuples = fromSeg2TGSeg(seg , ~0ul, false, 0, trackProvenance);
            const auto nodeId = (!shouldTrackProvenance() ||
                    tuples->getProvenanceType() == SEG_DIFFNODES) ? ~0ul :
                tuples->getNodeId();
            newtuples = outputJoin->getSegment(
                    nodeId,
                    false,
                    0,
                    tuples->getProvenanceType(),
                    extraColumns);
        }
    }
    return newtuples;
}

SegProvenanceType GBRuleExecutor::getSegProvenanceType() const {
    if (provenanceType == GBGraph::ProvenanceType::NOPROV) {
        return SegProvenanceType::SEG_NOPROV;
    } else if (provenanceType == GBGraph::ProvenanceType::NODEPROV) {
        return SegProvenanceType::SEG_SAMENODE;
    } else if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
        return SegProvenanceType::SEG_FULLPROV;
    } else {
        throw 10;
    }
}

std::vector<GBRuleOutput> GBRuleExecutor::executeRule(Rule &rule,
        GBRuleInput &node) {
    auto &bodyNodes = node.incomingEdges;

    lastDurationFirst = std::chrono::duration<double, std::milli>(0);
    lastDurationMergeSort = std::chrono::duration<double, std::milli>(0);
    lastDurationJoin = std::chrono::duration<double, std::milli>(0);
    lastDurationCreateHead = std::chrono::duration<double, std::milli>(0);
    bdyAtoms = "";

    //Perform the joins and populate the head
    auto &bodyAtoms = rule.getBody();
    //Maybe rearrange the body atoms? Don't forget to also re-arrange the body
    //nodes

#ifdef DEBUG
    //Check the cardinalities (this is not needed at the moment)
    size_t currentIDBBodyAtom = 0;
    std::vector<size_t> cardinalities;
    for(size_t i = 0; i < bodyAtoms.size(); ++i) {
        const Literal &currentBodyAtom = bodyAtoms[i];
        bool isEDB = currentBodyAtom.getPredicate().getType() == EDB;
        size_t c = 0;
        size_t nnodes = 1;
        if (isEDB) {
            c = layer.getCardinality(currentBodyAtom);
        } else {
            for(auto i : bodyNodes[currentIDBBodyAtom]) {
                c += g.getNodeData(i)->getNRows();
            }
            nnodes = bodyNodes[currentIDBBodyAtom].size();
            currentIDBBodyAtom++;
        }
        cardinalities.push_back(c);
        LOG(DEBUGL) << "Cardinality " << c << " nnodes " << nnodes;
    }

    currentIDBBodyAtom = 0;
    for(size_t i = 0; i < bodyAtoms.size(); ++i) {
        const Literal &currentBodyAtom = bodyAtoms[i];
        bool isEDB = currentBodyAtom.getPredicate().getType() == EDB;
        if (isEDB) {
            bdyAtoms += "-1 ";
        } else {
            auto nnodes = bodyNodes[currentIDBBodyAtom].size();
            bdyAtoms += std::to_string(nnodes) + " ";
            currentIDBBodyAtom++;
        }
    }
#endif

    std::vector<size_t> varsIntermediate;
    std::shared_ptr<const TGSegment> intermediateResults;
    //The following data structure is used only if trackProvenance=true
    std::vector<std::shared_ptr<Column>> intermediateResultsNodes;
    int64_t currentBodyNode = -1;
    int64_t prevBodyNode = -1;

    //Discover if some body atoms refer to built-in predicates that should not
    //be joined, but rather used as filtering functioning
    std::set<size_t> skippedBodyAtoms;
    std::vector<Literal> filteringBodyAtoms;
    for(size_t i = 0; i < bodyAtoms.size(); ++i) {
        const Literal &currentBodyAtom = bodyAtoms[i];
        bool isEDB = currentBodyAtom.getPredicate().getType() == EDB;
        if (isEDB) {
            if (!layer.acceptQueriesWithFreeVariables(currentBodyAtom)) {
                skippedBodyAtoms.insert(i);
                filteringBodyAtoms.push_back(currentBodyAtom);
            }
        }
    }

    bool enableCacheLeft = true;
    for(size_t i = 0; i < bodyAtoms.size(); ++i) {
        if (skippedBodyAtoms.count(i)) {
            enableCacheLeft = false;
            continue; //This atom is handled differently
        }
        const Literal &currentBodyAtom = bodyAtoms[i];
        VTuple currentVars = currentBodyAtom.getTuple();

        //Which are the positions (left,right) of the variable to join?
        std::vector<std::pair<int, int>> joinVarPos;
        //Which variables should be copied from the existing collection?
        std::vector<int> copyVarPosLeft;
        //Which variables should be copied from the new body atom?
        std::vector<int> copyVarPosRight;
        //Compute the value of the vars
        computeVarPos(varsIntermediate, i, bodyAtoms, rule.getHeads(),
                joinVarPos, copyVarPosLeft, copyVarPosRight);

        //Update the list of variables from the left atom
        std::vector<size_t> newVarsIntermediateResults;
        for(auto varIdx = 0; varIdx < copyVarPosLeft.size(); ++varIdx) {
            auto var = copyVarPosLeft[varIdx];
            newVarsIntermediateResults.push_back(varsIntermediate[var]);
        }
        //Update the list of variables from the right atom
        for(auto varIdx = 0; varIdx < copyVarPosRight.size(); ++varIdx) {
            auto varPos = copyVarPosRight[varIdx];
            newVarsIntermediateResults.push_back(
                    currentVars.get(varPos).getId());
        }

        //Compute possible builtin functions that would need to be applied
        //on the results of the join
        std::vector<BuiltinFunction> builtinFunctions;
        addBuiltinFunctions(builtinFunctions,
                filteringBodyAtoms, newVarsIntermediateResults);

        //Get the node of the current bodyAtom (if it's IDB)
        bool isCurrentBodyAtomEDB = currentBodyAtom.getPredicate().getType()
            == EDB;
        if (!isCurrentBodyAtomEDB) {
            currentBodyNode++;
        }
        if (i == 0) {
            //Process first body atom, no join
            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();
            if (isCurrentBodyAtomEDB) {
                enableCacheLeft = false;
                intermediateResults = processAtom_EDB(currentBodyAtom,
                        copyVarPosRight);
            } else {
                intermediateResults = processAtom_IDB(currentBodyAtom,
                        bodyNodes[currentBodyNode], copyVarPosRight, true, true);
            }
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            lastDurationFirst += end - start;
            durationFirst += lastDurationFirst;
            if (!builtinFunctions.empty()) {
                LOG(ERRORL) << "Builtin functions are not supported if they"
                    " use variables in one body atom";
                throw 10;
            }
            //If empty then stop
            if (intermediateResults == NULL || intermediateResults->isEmpty()) {
                intermediateResults = std::shared_ptr<const TGSegment>();
                break;
            }
        } else {
            std::vector<size_t> &nodesLeft = i == 1 && prevBodyNode >= 0 ?
                bodyNodes[prevBodyNode] : noBodyNodes;
            std::vector<size_t> &nodesRight =
                !isCurrentBodyAtomEDB ? bodyNodes[currentBodyNode] : noBodyNodes;

            uint8_t extraColumns = 0;
            bool multipleComb = provenanceType != GBGraph::ProvenanceType::NOPROV;
            if (multipleComb) {
                if (intermediateResults->getProvenanceType() == SEG_FULLPROV) {
                    size_t maxRightOffsetColumns = 1; //only offset, no node
                    extraColumns = 2 + intermediateResults->
                        getNOffsetColumns() - 1 + maxRightOffsetColumns;
                } else {
                    extraColumns = 2;
                }
            }
            bool removeDuplicates = retainUnique;
            std::unique_ptr<GBSegmentInserter> newIntermediateResults =
                GBSegmentInserter::getInserter(copyVarPosLeft.size() +
                        copyVarPosRight.size() + extraColumns,
                        extraColumns, removeDuplicates);
            newIntermediateResults->addBuiltinFunctions(builtinFunctions);

            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();

            bool enableCacheRight = true;
            if (currentBodyAtom.getNConstants() > 0 ||
                    !currentBodyAtom.getRepeatedVars().empty() ||
                    isCurrentBodyAtomEDB) {
                enableCacheRight = false;
            }

            //Perform a join (or a left join) between the intermediate results
            //and the new collection
            join(enableCacheLeft,
                    enableCacheRight,
                    intermediateResults,
                    nodesLeft,
                    nodesRight,
                    currentBodyAtom,
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    newIntermediateResults);
            enableCacheLeft = enableCacheLeft & enableCacheRight;

            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> durJoin = end - start;
            //LOG(INFOL) << "Join " << durJoin.count();
            lastDurationJoin += durJoin;
            durationJoin += durJoin;

            //If empty then stop
            if (newIntermediateResults->isEmpty()) {
                intermediateResults = std::shared_ptr<const TGSegment>();
                break;
            }

            if (shouldTrackProvenance()) {
                if  (multipleComb) {
                    //Process the output of nodes
                    newIntermediateResults->postprocessJoin(
                            intermediateResultsNodes, extraColumns);
                    intermediateResults = newIntermediateResults->getSegment(
                            ~0ul, false, 0, getSegProvenanceType(), extraColumns - 1);
                } else {
                    //This mode is activated only if there is at most one IDB
                    //atom in the body
                    size_t nodeId = ~0ul;
                    if (intermediateResults->getNodeId() != ~0ul) {
                        nodeId = intermediateResults->getNodeId();
                    }
                    if (!isCurrentBodyAtomEDB) {
                        assert(nodesRight.size() > 0);
                        nodeId = nodesRight[0];
                        assert(!g.isTmpNode(nodeId));
                    }
                    intermediateResults = newIntermediateResults->getSegment(
                            nodeId, false, 0, getSegProvenanceType(), 0);
                }
            } else {
                intermediateResults = newIntermediateResults->getSegment(
                        ~0ul, false, 0, getSegProvenanceType(), 0);
            }
        }
        if (!isCurrentBodyAtomEDB) {
            prevBodyNode++;
        }
        varsIntermediate = newVarsIntermediateResults;
    }

    //Filter out the derivations produced by the rule
    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    auto nonempty = !(intermediateResults == NULL ||
            intermediateResults->isEmpty());
    std::vector<GBRuleOutput> output;
    if (nonempty) {
        bool uniqueTuples = false;
        if (rule.isExistential()) {
            //Perform restricted check
            intermediateResults = performRestrictedCheck(rule,
                    intermediateResults, varsIntermediate);
            if (!intermediateResults->isEmpty()) {
                //If there are existential variables, add values for them
                intermediateResults = addExistentialVariables(rule,
                        intermediateResults, varsIntermediate);
            }
            uniqueTuples = true;
        }

        //Compute the head atoms
        if (!intermediateResults->isEmpty()) {
            for (auto &head : rule.getHeads()) {
                bool shouldSort = true, shouldDelDupl = true;
                shouldSortDelDupls(head, bodyAtoms, bodyNodes,
                        shouldSort, shouldDelDupl);
                auto results = projectHead(head,
                        varsIntermediate, intermediateResults,
                        shouldSort,
                        shouldDelDupl,
                        intermediateResultsNodes,
                        node);
                std::chrono::steady_clock::time_point end =
                    std::chrono::steady_clock::now();
                lastDurationCreateHead += end - start;
                durationCreateHead += lastDurationCreateHead;
                GBRuleOutput o;
                o.segment = results;
                o.nodes = intermediateResultsNodes;
                o.uniqueTuples = uniqueTuples;
                output.push_back(o);
            }
        }
    }
    return output;
}

void GBRuleExecutor::printStats() {
    LOG(INFOL) << "Time first (ms): " << durationFirst.count();
    LOG(INFOL) << "Time mergesort (ms): " << durationMergeSort.count();
    LOG(INFOL) << "Time joins (ms): " << durationJoin.count();
    LOG(INFOL) << "Time head (ms): " << durationCreateHead.count();
    //LOG(INFOL) << "Time preparation join2to1 (ms): " << durationPrep2to1.count();
}

DuplManager::DuplManager(GBSegmentInserter *output) {
    //Collect all the entities in the left side of the binary relation
    //added so far
    LOG(DEBUGL) << "Start populating the maps for the detection of duplicates";
    l = output->getEntitiesAddedSoFar(0);
    r = output->getEntitiesAddedSoFar(1);
    size_t nrows = output->getNRows();
    enabled = nrows >= (l.size() * r.size());
    LOG(DEBUGL) << "nrows = " << nrows << " l.size()=" <<  l.size() << " r.size()=" << r.size();
    LOG(DEBUGL) << "(Done) Start populating the maps for the detection of duplicates. Is enabled? " << enabled;
}

bool DuplManager::left(const Term_t &t) const {
    return enabled && l.count(t);
}

bool DuplManager::right(const Term_t &t) const {
    return enabled && r.count(t);
}
