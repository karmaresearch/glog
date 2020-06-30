#include <vlog/gbchase/gbruleexecutor.h>
#include <vlog/gbchase/gbsegmentcache.h>
#include <vlog/gbchase/gbsegment.h>
#include <vlog/gbchase/gbsegmentinserter.h>


std::shared_ptr<const TGSegment> GBRuleExecutor::projectTuples_structuresharing(
        std::shared_ptr<const TGSegment> tuples,
        const std::vector<int> &posKnownVariables,
        bool isSorted) {
    std::vector<std::shared_ptr<Column>> columns;
    tuples->projectTo(posKnownVariables, columns);
    return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(columns,
                tuples->getNRows(), isSorted, 0, trackProvenance));
}

std::chrono::duration<double, std::milli> GBRuleExecutor::getDuration(DurationType typ) {
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

std::shared_ptr<const TGSegment> GBRuleExecutor::projectTuples(
        std::shared_ptr<const TGSegment> tuples,
        const std::vector<int> &posKnownVariables) {
    if (posKnownVariables.size() == 0) {
        LOG(ERRORL) << "Projections with no columns are not supported";
        throw 10;
    } else if (posKnownVariables.size() == 1) {
        if (trackProvenance) {
            assert(tuples->getProvenanceType() != 0);
            if (tuples->getProvenanceType() == 1) {
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
        if (trackProvenance) {
            assert(tuples->getProvenanceType() != 0);
            if (tuples->getProvenanceType() == 1) {
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
    } else { //More than 2
        std::vector<std::shared_ptr<Column>> columns;
        tuples->projectTo(posKnownVariables, columns);
        return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(columns,
                    tuples->getNRows(), false, 0, trackProvenance));
    }
}

std::shared_ptr<const TGSegment> GBRuleExecutor::projectHead(
        const Literal &head,
        std::vector<size_t> &vars,
        std::shared_ptr<const TGSegment> tuples,
        bool shouldSort,
        bool shouldDelDupl) {
    const auto &tupleHead = head.getTuple();
    std::vector<int> posKnownVariables;
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
            LOG(ERRORL) << "Should not happen. At this stage, all variables "
                "are known";
            throw 10;
        }
    }

    //Project the variables in the order specified by the head atom
    assert(posKnownVariables.size() > 0); //For now, I do not support relations
    //with arity = 0
    if (posKnownVariables.size() == 1) {
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
            tuples = tuples->unique();
        }
    }
    return tuples;
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

std::shared_ptr<const TGSegment> GBRuleExecutor::processFirstAtom_EDB(
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

    bool shouldSort = true;
    bool shouldRetainUnique = true;

    auto &table = edbTables[p];
    std::vector<std::shared_ptr<Column>> columns;
    auto nrows = table->getCardinality(atom);
    if (table->useSegments()) {
        auto seg = table->getSegment();
        for(int i = 0; i < atom.getTupleSize(); ++i) {
            auto col = seg->getColumn(i);
            columns.push_back(col);
        }
    } else {
        //e.g., Trident database
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

    if (trackProvenance) {
        CompressedColumnBlock b(~0ul, 0, nrows);
        std::vector<CompressedColumnBlock> blocks;
        blocks.push_back(b);
        columns.push_back(std::shared_ptr<Column>(
                    new CompressedColumn(blocks, nrows)));
    }

    auto seg = std::shared_ptr<const TGSegment>(
            new TGSegmentLegacy(columns, nrows, true, 0, trackProvenance));
    if (copyVarPos.size() != atom.getTupleSize()) {
        //There is a projection. We might have to remove duplicates
        auto projectedSegment = projectTuples_structuresharing(
                seg, copyVarPos,
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
        return uniqueSegment;
    } else {
        return seg;
    }
}

void GBRuleExecutor::nestedloopjoin(
        std::shared_ptr<const TGSegment> inputLeft,
        const std::vector<size_t> &nodesLeft,
        std::shared_ptr<const TGSegment> inputRight,
        const std::vector<size_t> &nodesRight,
        const Literal &literalRight,
        std::vector<std::pair<int, int>> &joinVarsPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<GBSegmentInserter> &output) {

    std::vector<uint8_t> fields1;
    std::vector<uint8_t> fields2;
    for (uint32_t i = 0; i < joinVarsPos.size(); ++i) {
        fields1.push_back(joinVarsPos[i].first);
        fields2.push_back(joinVarsPos[i].second);
    }

    //Sort the left segment by the join variable
    if (!fields1.empty() && !inputLeft->isSortedBy(fields1)) {
        //Caching works only if there is at most one join variable...
        /*if (nodesLeft.size() > 0) {
          SegmentCache &c = SegmentCache::getInstance();
          if (!c.contains(nodesLeft, fields1)) {
          inputLeft = inputLeft->sortBy(fields1);
          c.insert(nodesLeft, fields1, inputLeft);
          } else {
          inputLeft = c.get(nodesLeft, fields1);
          }
          } else {*/
        inputLeft = inputLeft->sortBy(fields1);
        /*}*/
    }
    std::unique_ptr<TGSegmentItr> itrLeft = inputLeft->iterator();

    int64_t countLeft = 0;
    std::vector<Term_t> currentKey;
    VTuple t = literalRight.getTuple();
    std::vector<int> positions;
    for(int i = 0; i < t.getSize(); ++i) {
        positions.push_back(i);
    }
    auto sizerow = copyVarPosLeft.size() + copyVarPosRight.size();
    Term_t currentrow[sizerow + 2];

    while (itrLeft->hasNext()) {
        itrLeft->next();
        currentKey.clear();
        itrLeft->mark();
        for (int i = 0; i < fields1.size(); i++) {
            currentKey.push_back(itrLeft->get(fields1[i]));
        }
        //Do a lookup on the right-side
        for (int i = 0; i < fields2.size(); i++) {
            t.set(VTerm(0,currentKey[i]),fields2[i]);
        }
        Literal l(literalRight.getPredicate(), t);
        auto segRight = processFirstAtom_EDB(l, positions);
        auto itrRight = segRight->iterator();
        if (itrRight->hasNext()) {
            countLeft = 1; //First determine how many rows in the left side
            //share the same key
            while (itrLeft->hasNext()) {
                itrLeft->next();
                bool equal = true;
                for (int i = 0; i < fields1.size(); i++) {
                    auto k = itrLeft->get(fields1[i]);
                    if (k != currentKey[i]) {
                        equal = false;
                        break;
                    }
                }
                if (!equal) {
                    break;
                }
                countLeft++;
            }
            while (itrRight->hasNext()) {
                itrRight->next();
                size_t i = 0;
                itrLeft->reset();
                while (i < countLeft) {
                    //Materialize the join
                    for(int idx = 0; idx < copyVarPosLeft.size(); ++idx) {
                        auto leftPos = copyVarPosLeft[idx];
                        auto el = itrLeft->get(leftPos);
                        currentrow[idx] = el;
                    }
                    for(int idx = 0; idx < copyVarPosRight.size(); ++idx) {
                        auto rightPos = copyVarPosRight[idx];
                        auto value = itrRight->get(rightPos);
                        currentrow[copyVarPosLeft.size() + idx] = value;
                    }
                    if (trackProvenance) {
                        currentrow[sizerow] = itrLeft->getNodeId();
                        currentrow[sizerow + 1] = itrRight->getNodeId();
                    }
                    output->add(currentrow);
                    if (i < (countLeft - 1))
                        itrLeft->next();
                    i++;
                }
            }

        }
    }
}

void GBRuleExecutor::joinTwoOne_EDB(
        std::shared_ptr<const TGSegment> inputLeft,
        std::shared_ptr<const TGSegment> inputRight,
        int joinLeftVarPos,
        std::vector<int> &copyVarPosLeft,
        std::unique_ptr<GBSegmentInserter> &output) {
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();

    //Get column 1, 2 left
    auto col1Left = ((TGSegmentLegacy*)inputLeft.get())
        ->getColumn(0);
    auto col2Left = ((TGSegmentLegacy*)inputLeft.get())
        ->getColumn(1);

    //Get column 1 right
    auto col1Right = ((TGSegmentLegacy*)inputRight.get())
        ->getColumn(0);

    if (!col1Left->isEDB() || !col2Left->isEDB() || !col1Right->isEDB()) {
        LOG(ERRORL) << "Not supported";
        throw 10;
    }

    const Literal &l1Left = ((EDBColumn*)col1Left.get())->getLiteral();
    uint8_t pos1 = ((EDBColumn*)col1Left.get())->posColumnInLiteral();
    std::vector<uint8_t> posColumns1;
    posColumns1.push_back(pos1);

    const Literal &l2Left = ((EDBColumn*)col2Left.get())->getLiteral();
    pos1 = ((EDBColumn*)col2Left.get())->posColumnInLiteral();
    posColumns1.push_back(pos1);

    std::vector<Substitution> subs;
    if (!l1Left.sameVarSequenceAs(l2Left) ||
            l2Left.subsumes(subs, l1Left, l2Left) == -1) {
        //The columns come from different literals. This is not yet
        //supported
        LOG(ERRORL) << "Not supported";
        throw 10;
    }

    const Literal &l1Right = ((EDBColumn*)col1Right.get())->getLiteral();
    uint8_t pos2 = ((EDBColumn*)col1Right.get())->posColumnInLiteral();

    const uint8_t ncopyvars = copyVarPosLeft.size();
    if (ncopyvars == 1) {
        std::vector<Term_t> c1;
        layer.join(c1,
                l1Left, posColumns1, joinLeftVarPos, l1Right, pos2,
                copyVarPosLeft[0]);
        output->swap(c1);
    } else if (ncopyvars == 2) {
        std::vector<std::pair<Term_t,Term_t>> c2;
        layer.join(c2,
                l1Left, posColumns1, joinLeftVarPos, l1Right, pos2,
                copyVarPosLeft[0], copyVarPosLeft[1]);
        output->swap(c2);
    } else {
        LOG(ERRORL) << "Case not supported";
        throw 10;
    }

    //dur = std::chrono::system_clock::now() - start;
    //LOG(INFOL) << "2To1 TotalJoin: " << dur.count();
    //LOG(INFOL) << "2To1 Copy: " << durationCopy.count();
}

void GBRuleExecutor::joinTwoOne(
        std::shared_ptr<const TGSegment> inputLeft,
        std::shared_ptr<const TGSegment> inputRight,
        int joinLeftVarPos,
        std::vector<int> &copyVarPosLeft,
        const int copyNodes,
        std::unique_ptr<GBSegmentInserter> &output) {
    //Preparation
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
    std::vector<uint8_t> fields1;
    fields1.push_back(joinLeftVarPos);
    std::vector<uint8_t> fields2;
    fields2.push_back(0);
    inputLeft = inputLeft->sortBy(fields1);
    std::unique_ptr<TGSegmentItr> itrLeft = inputLeft->iterator();
    inputRight = inputRight->sortBy(fields2);
    std::unique_ptr<TGSegmentItr> itrRight = inputRight->iterator();
    std::chrono::duration<double, std::milli> dur =
        std::chrono::system_clock::now() - start;
    //lastDurationPrep2to1 += dur;
    //durationPrep2to1 += dur;

    const uint8_t ncopyvars = copyVarPosLeft.size();
    const uint8_t varpos1 = copyVarPosLeft[0];
    const uint8_t varpos2 = ncopyvars == 2 ? copyVarPosLeft[1] : 0;

    const bool fastCopy = copyNodes == 0 &&
        output->getNBuiltinFunctions() == 0;
    std::vector<Term_t> c1;
    std::vector<std::pair<Term_t,Term_t>> c2;
    Term_t currentrow[4];

    std::chrono::duration<double, std::milli> durationCopy;

    start = std::chrono::system_clock::now();
    //Join
    Term_t v1 = ~0ul;
    Term_t v2 = ~0ul;
    while (true) {
        if (v1 == ~0ul) {
            if (itrLeft->hasNext()) {
                itrLeft->next();
                v1 = itrLeft->get(joinLeftVarPos);
            } else {
                break;
            }
        }
        if (v2 == ~0ul) {
            if (itrRight->hasNext()) {
                itrRight->next();
                v2 = itrRight->get(0);
            } else {
                break;
            }
        }
        if (v1 == v2) {
            std::chrono::system_clock::time_point start1 =
                std::chrono::system_clock::now();
            if (fastCopy) {
                if (ncopyvars == 1) {
                    c1.push_back(itrLeft->get(varpos1));
                } else {
                    c2.push_back(std::make_pair(itrLeft->get(varpos1),
                                itrLeft->get(varpos2)));
                }
            } else {
                currentrow[0] = itrLeft->get(varpos1);
                int startidx = 1;
                if (ncopyvars == 2) {
                    currentrow[1] = itrLeft->get(varpos2);
                    startidx = 2;
                }
                if (copyNodes == 2) {
                    currentrow[startidx] = itrRight->getNodeId();
                    currentrow[startidx + 1] = itrLeft->getNodeId();
                } else if (copyNodes == 1) {
                    currentrow[startidx] = itrLeft->getNodeId();
                    currentrow[startidx + 1] = itrRight->getNodeId();
                }
                output->add(currentrow);
            }
            v1 = ~0ul;
            auto dur1 = std::chrono::system_clock::now() - start1;
            durationCopy += dur1;
        } else if (v1 < v2) {
            v1 = ~0ul;
        } else {
            v2 = ~0ul;
        }
    }

    if (fastCopy && ncopyvars == 1) {
        output->swap(c1);
    }
    if (fastCopy && ncopyvars == 2) {
        output->swap(c2);
    }

    //dur = std::chrono::system_clock::now() - start;
    //LOG(INFOL) << "2To1 TotalJoin: " << dur.count();
    //LOG(INFOL) << "2To1 Copy: " << durationCopy.count();
}

void GBRuleExecutor::mergejoin(
        std::shared_ptr<const TGSegment> inputLeft,
        const std::vector<size_t> &nodesLeft,
        std::shared_ptr<const TGSegment> inputRight,
        const std::vector<size_t> &nodesRight,
        std::vector<std::pair<int, int>> &joinVarsPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<GBSegmentInserter> &output) {
    std::chrono::system_clock::time_point startL =
        std::chrono::system_clock::now();

    std::vector<uint8_t> fields1;
    std::vector<uint8_t> fields2;
    for (uint32_t i = 0; i < joinVarsPos.size(); ++i) {
        fields1.push_back(joinVarsPos[i].first);
        fields2.push_back(joinVarsPos[i].second);
    }
    std::chrono::system_clock::time_point startS =
        std::chrono::system_clock::now();

    //Sort the left segment by the join variable
    if (!fields1.empty() && !inputLeft->isSortedBy(fields1)) {
        if (fields1.size() > 1)
            LOG(WARNL) << "I cannot use the cache because it works only if the join is over 1 variable";
        if (nodesLeft.size() > 0 && fields1.size() == 1) {
            SegmentCache &c = SegmentCache::getInstance();
            if (!c.contains(nodesLeft, fields1)) {
                inputLeft = inputLeft->sortBy(fields1);
                c.insert(nodesLeft, fields1, inputLeft);
            } else {
                inputLeft = c.get(nodesLeft, fields1);
            }
        } else {
            inputLeft = inputLeft->sortBy(fields1);
        }
    }
    std::unique_ptr<TGSegmentItr> itrLeft = inputLeft->iterator();
    //Sort the right segment by the join variable
    if (!fields2.empty() && !inputRight->isSortedBy(fields2)) {
        if (fields2.size() > 1)
            LOG(WARNL) << "I cannot use the cache because it works only if the join is over 1 variable";
        if (nodesRight.size() > 0 && fields2.size() == 1) {
            SegmentCache &c = SegmentCache::getInstance();
            if (fields2.size() == 1) {
                if (!c.contains(nodesRight, fields2)) {
                    inputRight = inputRight->sortBy(fields2);
                    c.insert(nodesRight, fields2, inputRight);
                } else {
                    inputRight = c.get(nodesRight, fields2);
                }
            } else {
                assert(fields2.size() > 1);
                inputRight = inputRight->sortBy(fields2);
            }
        } else {
            inputRight = inputRight->sortBy(fields2);
        }
    }
    std::unique_ptr<TGSegmentItr> itrRight = inputRight->iterator();
    lastDurationMergeSort += std::chrono::system_clock::now() - startS;
    durationMergeSort += lastDurationMergeSort;

#if NDEBUG
    size_t joinLeftSize = inputLeft->getNRows();
    size_t joinRightSize = inputRight->getNRows();
    LOG(DEBUGL) << "Join left size=" << joinLeftSize << " right size=" << joinRightSize;
#endif

    //Do the merge join
    if (itrLeft->hasNext()) {
        itrLeft->next();
    } else {
        return;
    }

    if (itrRight->hasNext()) {
        itrRight->next();
    } else {
        return;
    }

#if NDEBUG
    size_t total = 0;
    size_t max = 65536;
#endif

    long countLeft = -1;
    std::vector<Term_t> currentKey;
    auto sizerow = copyVarPosLeft.size() + copyVarPosRight.size();
    Term_t currentrow[sizerow + 2];
    for(int i = 0; i < sizerow + 2; ++i) currentrow[i] = 0;
    int res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(), joinVarsPos);
    size_t processedRight = 0;
    while (true) {
        //Are they matching?
        while (res < 0 && itrLeft->hasNext()) {
            itrLeft->next();
            res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(), joinVarsPos);
        }

#if NDEBUG
        if (processedRight % 1000 == 0)
            LOG(DEBUGL) << "Processed records " << processedRight;
#endif

        if (res < 0) //The first iterator is finished
            break;

        while (res > 0 && itrRight->hasNext()) {
            itrRight->next();
#if NDEBUG
            processedRight++;
#endif
            res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(), joinVarsPos);
        }

        if (res > 0) { //The second iterator is finished
            break;
        }

        if (res == 0) {
            if (countLeft == -1) {
                currentKey.clear();
                itrLeft->mark();
                countLeft = 1;
                for (int i = 0; i < fields1.size(); i++) {
                    currentKey.push_back(itrLeft->get(fields1[i]));
                }
                while (itrLeft->hasNext()) {
                    itrLeft->next();
                    bool equal = true;
                    for (int i = 0; i < fields1.size(); i++) {
                        auto k = itrLeft->get(fields1[i]);
                        if (k != currentKey[i]) {
                            equal = false;
                            break;
                        }
                    }
                    if (! equal) {
                        break;
                    }
                    countLeft++;
                }
            }
#if NDEBUG
            total += countLeft;
            while (total >= max) {
                LOG(TRACEL) << "Count = " << countLeft << ", total = " << total;
                max = max + max;
            }
#endif

            itrLeft->reset();
            //Move the left iterator countLeft times and emit tuples
            for(int idx = 0; idx < copyVarPosRight.size(); ++idx) {
                auto rightPos = copyVarPosRight[idx];
                currentrow[copyVarPosLeft.size() + idx] = itrRight->get(
                        rightPos);
            }
            auto c = 0;
            bool leftActive = true;
            while (c < countLeft) {
                for(int idx = 0; idx < copyVarPosLeft.size(); ++idx) {
                    auto leftPos = copyVarPosLeft[idx];
                    auto el = itrLeft->get(leftPos);
                    currentrow[idx] = el;
                }
                if (trackProvenance) {
                    currentrow[sizerow] = itrLeft->getNodeId();
                    currentrow[sizerow + 1] = itrRight->getNodeId();
                }
                output->add(currentrow);
                if (itrLeft->hasNext()) {
                    itrLeft->next();
                } else {
                    leftActive = false;
                }
                c++;
            }

            //Move right
            if (!itrRight->hasNext()) {
                break;
            } else {
#if NDEBUG
                processedRight++;
#endif
                itrRight->next();
            }
            //Does the right row have not the same value as before?
            bool equal = true;
            for (int i = 0; i < fields2.size(); i++) {
                auto newKey = itrRight->get(fields2[i]);
                if (newKey != currentKey[i]) {
                    equal = false;
                    break;
                }
            }
            if (! equal) {
                countLeft = -1;
                //LeftItr is already pointing to the next element ...
                if (!leftActive) {
                    break;
                }
                res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(),
                        joinVarsPos);
            }
        }
    }
#if NDEBUG
    LOG(DEBUGL) << "Total = " << total;
    std::chrono::duration<double> secL =
        std::chrono::system_clock::now() - startS;
    LOG(TRACEL) << "merge_join: time : " << secL.count() * 1000;
#endif
}

void GBRuleExecutor::leftjoin(
        std::shared_ptr<const TGSegment> inputLeft,
        const std::vector<size_t> &nodesLeft,
        std::shared_ptr<const TGSegment> inputRight,
        std::vector<std::pair<int, int>> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::unique_ptr<GBSegmentInserter> &output,
        const bool copyOnlyLeftNode) {

    std::vector<uint8_t> fields1;
    std::vector<uint8_t> fields2;
    for (uint32_t i = 0; i < joinVarPos.size(); ++i) {
        fields1.push_back(joinVarPos[i].first);
        fields2.push_back(joinVarPos[i].second);
    }

    //Sort the left segment by the join variable
    if (!fields1.empty() && !inputLeft->isSortedBy(fields1)) {
        if (nodesLeft.size() > 0) {
            SegmentCache &c = SegmentCache::getInstance();
            if (!c.contains(nodesLeft, fields1)) {
                inputLeft = inputLeft->sortBy(fields1);
                c.insert(nodesLeft, fields1, inputLeft);
            } else {
                inputLeft = c.get(nodesLeft, fields1);
            }
        } else {
            inputLeft = inputLeft->sortBy(fields1);
        }
    }
    auto itrLeft = inputLeft->iterator();

    //Get all the triples in the right relation
    auto sortedInputRight = inputRight->sortBy(fields2);
    auto itrRight = sortedInputRight->iterator();

    bool leftActive = false;
    bool rightActive = false;
    if (itrLeft->hasNext()) {
        itrLeft->next();
        leftActive = true;
    } else {
        return;
    }
    if (itrRight->hasNext()) {
        itrRight->next();
        rightActive = true;
    }

    auto sizerow = copyVarPosLeft.size();
    Term_t currentrow[sizerow + 2];

    while (leftActive && rightActive) {
        int res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(), joinVarPos);
        if (res <= 0) {
            if (res < 0) {
                for(int idx = 0; idx < copyVarPosLeft.size(); ++idx) {
                    auto leftPos = copyVarPosLeft[idx];
                    auto el = itrLeft->get(leftPos);
                    currentrow[idx] = el;
                }
                output->add(currentrow);
            }
            if (itrLeft->hasNext()) {
                itrLeft->next();
            } else {
                leftActive = false;
            }
        } else {
            if (itrRight->hasNext()) {
                itrRight->next();
            } else {
                rightActive = false;
            }
        }
    }
    while (leftActive) {
        for(int idx = 0; idx < copyVarPosLeft.size(); ++idx) {
            auto leftPos = copyVarPosLeft[idx];
            auto el = itrLeft->get(leftPos);
            currentrow[idx] = el;
        }
        if (trackProvenance) {
            currentrow[sizerow] = itrLeft->getNodeId();
            if (!copyOnlyLeftNode)
                currentrow[sizerow + 1] = 0;
        }
        output->add(currentrow);
        if (itrLeft->hasNext()) {
            itrLeft->next();
        } else {
            leftActive = false;
        }
    }
}

void GBRuleExecutor::join(
        std::shared_ptr<const TGSegment> inputLeft,
        const std::vector<size_t> &nodesLeft,
        std::vector<size_t> &nodesRight,
        const Literal &literalRight,
        std::vector<std::pair<int, int>> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<GBSegmentInserter> &output) {

    std::shared_ptr<const TGSegment> inputRight;
    bool mergeJoinPossible = true;
    if (nodesRight.size() == 1) {
        size_t idbBodyAtomIdx = nodesRight[0];
        inputRight = g.getNodeData(idbBodyAtomIdx);
    } else if (nodesRight.size() > 0) {
        auto ncols = g.getNodeData(nodesRight[0])->getNColumns();
        std::vector<int> projectedPos;
        for(int i = 0; i < ncols; ++i)
            projectedPos.push_back(i);
        inputRight = g.mergeNodes(nodesRight, projectedPos);
    } else {
        //It must be an EDB literal because only these do not have nodes
        assert(literalRight.getPredicate().getType() == EDB);

        //check if I'm allowed to do a merge join
        if (!layer.isQueryAllowed(literalRight)) {
            mergeJoinPossible = false;
        } else {
            std::vector<int> allVars;
            for(int i = 0; i < literalRight.getTupleSize(); ++i)
                allVars.push_back(i);
            inputRight = processFirstAtom_EDB(literalRight, allVars);
        }
    }

    if (literalRight.isNegated()) {
        // Negated atoms should not introduce new variables.
        assert(copyVarPosRight.size() == 0);
        leftjoin(inputLeft,
                nodesLeft,
                inputRight,
                joinVarPos,
                copyVarPosLeft,
                output);
    } else {
        if (mergeJoinPossible) {
            /*if (inputLeft->getNColumns() == 2
              && inputRight->getNColumns() == 1
              && joinVarPos.size() == 1) {
              assert(copyVarPosRight.size() == 0);
              int copyTypeNode = 0;
              if (inputLeft->getProvenanceType() == 2 ||
              inputRight->getProvenanceType() == 2) {
              copyTypeNode = 1;
              }
              const uint8_t jv = joinVarPos[0].first;

              if (!copyTypeNode && inputLeft->hasColumnarBackend()
              && inputRight->hasColumnarBackend()) {
            //I can speed up the join with a join directly at the EDB
            //level
            joinTwoOne_EDB(inputLeft,
            inputRight,
            jv,
            copyVarPosLeft,
            output);
            } else {
            joinTwoOne(
            inputLeft,
            inputRight,
            jv,
            copyVarPosLeft,
            copyTypeNode,
            output);
            }
            } else if (inputLeft->getNColumns() == 1
            && inputRight->getNColumns() == 2
            && joinVarPos.size() == 1) {
            int copyTypeNode = 0;
            if (inputLeft->getProvenanceType() == 2 ||
            inputRight->getProvenanceType() == 2)
            copyTypeNode = 2;
            const uint8_t jv = joinVarPos[0].second;
            std::vector<int> cp;
            if (copyVarPosLeft.size() == 1) {
            cp.push_back(jv);
            }
            assert(copyVarPosRight.size() < 2);
            if (copyVarPosRight.size() > 0) {
            cp.push_back(copyVarPosRight[0]);
            }

            if (!copyTypeNode && inputLeft->hasColumnarBackend()
            && inputRight->hasColumnarBackend()) {
            //I can speed up the join with a join directly at the EDB
            //level
            joinTwoOne_EDB(
            inputRight,
            inputLeft,
            jv,
            cp,
            output);
            } else {
            joinTwoOne(
            inputRight,
            inputLeft,
            jv,
            cp,
            copyTypeNode,
            output);
            }
            } else {*/
            mergejoin(inputLeft,
                    nodesLeft,
                    inputRight,
                    nodesRight,
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    output);
            //}
        } else {
            nestedloopjoin(inputLeft,
                    nodesLeft,
                    inputRight,
                    nodesRight,
                    literalRight,
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    output);
        }
    }
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
        bool sortOk = true;
        for(int i = 0; i < th.getSize(); ++i) {
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
                }
            }
        }
        shouldSort = (!sortOk) ||
            (bodyNodes.size() > 0 && bodyNodes[0].size() > 1);


    } else {
        shouldSort = shouldDelDupl = true;
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
            mode = !trackProvenance ? 0 : tuples->getProvenanceType() == 1 ?
                1 : 2;
        } else if (nfields == 2) {
            mode = !trackProvenance ? 3 : tuples->getProvenanceType() == 1 ?
                4 : 5;
        } else {
            if (!trackProvenance) {
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
                        new UnaryWithConstProvTGSegment(terms1, ~0ul, false,
                            0));
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
                        new BinaryWithConstProvTGSegment(terms2, ~0ul, false,
                            0));
                break;
            case 5:
                tuples = std::shared_ptr<const TGSegment>(
                        new BinaryWithProvTGSegment(terms3, ~0ul, false, 0));
                break;
            case 6:
            case 7:
                auto seg = terms4->getSegment();
                std::vector<std::shared_ptr<Column>> columns;
                auto nfieldsToCopy = (trackProvenance) ? nfields + 1 : nfields;
                for(int i = 0; i < nfieldsToCopy; ++i) {
                    columns.push_back(seg->getColumn(i));
                }
                tuples = std::shared_ptr<const TGSegment>(
                        new TGSegmentLegacy(columns, seg->getNRows(), false,
                            0, mode == 7));
                break;
        }
    }
    g.setCounterNullValues(counter);
    return tuples;
}

std::shared_ptr<const TGSegment> GBRuleExecutor::performRestrictedCheck(
        Rule &rule,
        std::shared_ptr<const TGSegment> tuples,
        const std::vector<size_t> &varTuples) {
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
            const int nfields = trackProvenance ? tuples->getNColumns() + 1 :
                tuples->getNColumns();
            std::unique_ptr<GBSegmentInserter> outputJoin = GBSegmentInserter::
                getInserter(nfields);

            //Perform a left join
            leftjoin(tuples,
                    nodes,
                    inputRight,
                    joinVarPos,
                    copyVarPosLeft,
                    outputJoin,
                    true);

            //Create a TGSegment from SegmentInserter
            //std::shared_ptr<const Segment> seg = outputJoin->getSegment();
            //tuples = fromSeg2TGSeg(seg , ~0ul, false, 0, trackProvenance);
            tuples = outputJoin->getSegment(~0ul, false, 0, trackProvenance);
        }
    }
    return tuples;
}

std::vector<GBRuleOutput> GBRuleExecutor::executeRule(Rule &rule,
        GBRuleInput &node) {
    auto &bodyNodes = node.incomingEdges;

#ifdef DEBUG//DEBUG
    if (rule.getFrontierVariables().empty()) {
        LOG(ERRORL) << "The system does not yet support the execution of rules"
            " with empty frontier variables set";
        throw 10;
    }
    LOG(DEBUGL) << "Execute rule " << node.ruleIdx << " " << rule.tostring(program, &layer);
#endif

    lastDurationFirst = std::chrono::duration<double, std::milli>(0);
    lastDurationMergeSort = std::chrono::duration<double, std::milli>(0);
    lastDurationJoin = std::chrono::duration<double, std::milli>(0);
    lastDurationCreateHead = std::chrono::duration<double, std::milli>(0);
    //lastDurationPrep2to1 = std::chrono::duration<double, std::milli>(0);
    //bdyAtoms = "";

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
        LOG(INFOL) << "Cardinality " << c << " nnodes " << nnodes;
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
    size_t currentBodyNode = 0;
    bool firstBodyAtomIsIDB = false;

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

    for(size_t i = 0; i < bodyAtoms.size(); ++i) {
        if (skippedBodyAtoms.count(i)) {
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
        bool isEDB = currentBodyAtom.getPredicate().getType() == EDB;

        if (i == 0) {
            //Process first body atom, no join
            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();
            if (isEDB) {
                intermediateResults = processFirstAtom_EDB(currentBodyAtom,
                        copyVarPosRight);
            } else {
                //if (bodyNodes[currentBodyNode].size() == 1) {
                //    auto n = bodyNodes[currentBodyNode][0];
                //}
                firstBodyAtomIsIDB = true;
                intermediateResults = g.mergeNodes(
                        bodyNodes[currentBodyNode], copyVarPosRight);
                currentBodyNode++;
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
            if (intermediateResults->isEmpty()) {
                intermediateResults = std::shared_ptr<const TGSegment>();
                break;
            }
        } else {
            std::vector<size_t> &nodesRight =
                !isEDB ? bodyNodes[currentBodyNode] : noBodyNodes;
            uint8_t extraColumns = 0;
            bool multipleComb = intermediateResults->getProvenanceType() == 2 ||
                nodesRight.size() > 1;
            if (trackProvenance && multipleComb) {
                extraColumns = 2;
            }
            std::unique_ptr<GBSegmentInserter> newIntermediateResults =
                GBSegmentInserter::getInserter(copyVarPosLeft.size() +
                        copyVarPosRight.size() + extraColumns);
            newIntermediateResults->addBuiltinFunctions(builtinFunctions);

            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();

            //Perform a join (or a left join) between the intermediate results
            //and the new collection
            join(intermediateResults,
                    i == 1 && firstBodyAtomIsIDB ? bodyNodes[0] : noBodyNodes,
                    nodesRight,
                    currentBodyAtom,
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    newIntermediateResults);

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

            if (trackProvenance) {
                if (multipleComb) {
                    //Process the output of nodes
                    newIntermediateResults->postprocessJoin(
                            intermediateResultsNodes);
                    intermediateResults = newIntermediateResults->getSegment(
                            ~0ul, false, 0, trackProvenance, true);
                } else {
                    //Add two columns in intermediateResultsNodes
                    intermediateResultsNodes.push_back(std::shared_ptr<Column>(
                                new CompressedColumn(
                                    intermediateResults->getNodeId(),
                                    newIntermediateResults->getNRows())));
                    size_t secondNodeId = ~0ul;
                    if (nodesRight.size() == 1) {
                        secondNodeId = g.getNodeData(nodesRight[0])->getNodeId();
                    }
                    assert(!g.isTmpNode(secondNodeId));
                    intermediateResultsNodes.push_back(std::shared_ptr<Column>(
                                new CompressedColumn(
                                    secondNodeId,
                                    newIntermediateResults->getNRows())));
                    intermediateResults = newIntermediateResults->getSegment(
                            ~0ul, false, 0, trackProvenance, false);
                }
            } else {
                intermediateResults = newIntermediateResults->getSegment(
                        ~0ul, false, 0, trackProvenance, false);
            }
            currentBodyNode++;
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
            //If there are existential variables, add values for them
            intermediateResults = addExistentialVariables(rule,
                    intermediateResults, varsIntermediate);
            uniqueTuples = true;
        }

        //If there was only one body atom, then the vector with the provenance
        //is empty
        /*        if (trackProvenance && firstBodyAtomIsIDB &&
                  bodyAtoms.size() == 1) {
                  size_t ni = bodyNodes[0][0];
                  if (g.isTmpNode(ni))
                  ni = intermediateResults->getNodeId();
        //assert(ni == intermediateResults->getNodeId());
        intermediateResultsNodes.push_back(std::shared_ptr<Column>(
        new CompressedColumn(ni,
        intermediateResults->getNRows())));
        }*/

        //Compute the head atoms
        for (auto &head : rule.getHeads()) {
            bool shouldSort = true, shouldDelDupl = true;
            shouldSortDelDupls(head, bodyAtoms, bodyNodes,
                    shouldSort, shouldDelDupl);
            auto results = projectHead(head,
                    varsIntermediate, intermediateResults,
                    shouldSort,
                    shouldDelDupl);
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
    return output;
}

void GBRuleExecutor::printStats() {
    LOG(INFOL) << "Time first (ms): " << durationFirst.count();
    LOG(INFOL) << "Time mergesort (ms): " << durationMergeSort.count();
    LOG(INFOL) << "Time joins (ms): " << durationJoin.count();
    LOG(INFOL) << "Time head (ms): " << durationCreateHead.count();
    //LOG(INFOL) << "Time preparation join2to1 (ms): " << durationPrep2to1.count();
}
