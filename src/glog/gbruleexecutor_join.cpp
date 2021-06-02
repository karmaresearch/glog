#include <glog/gbruleexecutor.h>
#include <glog/gbsegmentcache.h>

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
    if (nodesRight.size() == 1 && provenanceType != GBGraph::ProvenanceType::FULLPROV) {
        size_t idbBodyAtomIdx = nodesRight[0];
        inputRight = g.getNodeData(idbBodyAtomIdx);
    } else if (nodesRight.size() > 0) {
        auto ncols = g.getNodeData(nodesRight[0])->getNColumns();
        std::vector<int> projectedPos;
        for(int i = 0; i < ncols; ++i)
            projectedPos.push_back(i);
        if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
            inputRight = processAtom_IDB(literalRight, nodesRight, projectedPos, false, true);
        } else {
            inputRight = processAtom_IDB(literalRight, nodesRight, projectedPos, false, false);
        }
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
            inputRight = processAtom_EDB(literalRight, allVars);
        }
    }
    if (inputRight->getNRows() == 0)
        return;

    if (literalRight.isNegated()) {
        if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
            LOG(ERRORL) << "Left joins not supported in FULLPROV mode";
            throw 10;
        }

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
            mergejoin(inputLeft,
                    nodesLeft,
                    inputRight,
                    nodesRight,
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    output);
        } else {
            if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
                LOG(ERRORL) << "Nested joins not supported in FULLPROV mode";
                throw 10;
            }
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
            LOG(WARNL) << "I cannot use the cache because it works only if "
                "the join is over 1 variable";
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
            LOG(WARNL) << "I cannot use the cache because it works only if "
                "the join is over 1 variable";
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

#if DEBUG
    size_t joinLeftSize = inputLeft->getNRows();
    size_t joinRightSize = inputRight->getNRows();
    LOG(DEBUGL) << "Join left size=" << joinLeftSize << " right size=" << joinRightSize;
#endif

    //START Datastructures used in case the join produces many duplicates
    size_t countDuplicatedJoins = 0;
    size_t addedSoFar = 0;
    size_t grpCountRight = 0;
    bool filterDuplEnabled = false;
    std::unique_ptr<DuplManager> duplManager;
    bool leftGroupDuplicate = false;
    bool rightGroupDuplicate = false;
    //END data structures used with many duplicates

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

#if DEBUG
    size_t total = 0;
    size_t max = 65536;
    size_t processedRight = 0;
#endif

    long countLeft = -1;
    std::vector<Term_t> currentKey;
    auto extraLeft = 0;
    if (inputLeft->getNOffsetColumns() > 0)
        extraLeft = inputLeft->getNOffsetColumns() - 1; //I remove the node,
    //it will be added at the end
    auto sizeLeftSide = copyVarPosLeft.size() + extraLeft;
    auto extraRight = 0;
    if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
        extraRight = 1;
    }
    auto sizeRightSide = copyVarPosRight.size() + extraRight;

    auto sizeRow = sizeLeftSide + sizeRightSide;
    Term_t currentrow[sizeRow + 2];
    for(size_t i = 0; i < sizeRow + 2; ++i) currentrow[i] = 0;
    int res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(), joinVarsPos);
    while (true) {
        //Are they matching?
        while (res < 0 && itrLeft->hasNext()) {
            itrLeft->next();
            res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(), joinVarsPos);
        }

#if DEBUG
        if (processedRight % 100000 == 0)
            LOG(DEBUGL) << "Processed records " << processedRight;
#endif

        if (res < 0) //The first iterator is finished
            break;

        while (res > 0 && itrRight->hasNext()) {
            itrRight->next();
#if DEBUG
            processedRight++;
#endif
            res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(), joinVarsPos);
        }

        if (res > 0) { //The second iterator is finished
            break;
        }

        if (res == 0) {
            if (countLeft == -1) {
                //ignored if dupl. are not considered
                leftGroupDuplicate = filterDuplEnabled;

                currentKey.clear();
                itrLeft->mark();
                countLeft = 1;
                grpCountRight = 1;

                for (int i = 0; i < fields1.size(); i++) {
                    auto v = itrLeft->get(fields1[i]);
                    currentKey.push_back(v);
                    leftGroupDuplicate = filterDuplEnabled &&
                        duplManager->left(itrLeft->get(copyVarPosLeft[0]));
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
                    if (!equal) {
                        break;
                    }
                    //Check for duplicates
                    if (leftGroupDuplicate) {
                        leftGroupDuplicate = filterDuplEnabled &&
                            duplManager->left(itrLeft->get(copyVarPosLeft[0]));
                    }
                    countLeft++;
                }
            }
#if DEBUG
            total += countLeft;
            while (total >= max) {
                LOG(TRACEL) << "Count = " << countLeft << ", total = " << total;
                max = max + max;
            }
#endif

            bool leftActive = true;
            rightGroupDuplicate = filterDuplEnabled &&
                duplManager->right(itrRight->get(copyVarPosRight[0]));
            if (!leftGroupDuplicate || !rightGroupDuplicate) {
                //Move the left iterator countLeft times and emit tuples
                for(int idx = 0; idx < copyVarPosRight.size(); ++idx) {
                    auto rightPos = copyVarPosRight[idx];
                    currentrow[copyVarPosLeft.size() + idx] = itrRight->get(rightPos);
                }
                for(int idx = 0; idx < extraRight; ++idx) {
                    currentrow[copyVarPosLeft.size() + copyVarPosRight.size()
                        + extraLeft + idx] = itrRight->getProvenanceOffset(idx);
                }
                itrLeft->reset();
                size_t c = 0;
                while (c < countLeft) {
                    for(int idx = 0; idx < copyVarPosLeft.size(); ++idx) {
                        auto leftPos = copyVarPosLeft[idx];
                        auto el = itrLeft->get(leftPos);
                        currentrow[idx] = el;
                    }
                    for(int idx = 0; idx < extraLeft; ++idx) {
                        currentrow[copyVarPosLeft.size() + copyVarPosRight.size() + idx] =
                            itrLeft->getProvenanceOffset(idx);
                    }
                    if (shouldTrackProvenance()) {
                        currentrow[sizeRow] = itrLeft->getNodeId();
                        currentrow[sizeRow + 1] = itrRight->getNodeId();
                    }
                    output->add(currentrow);
                    if (itrLeft->hasNext()) {
                        itrLeft->next();
                    } else {
                        leftActive = false;
                    }
                    c++;
                }
            }

            //Move right
            if (!itrRight->hasNext()) {
                break;
            } else {
#if DEBUG
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
            if (!equal) {
                //The algorithm should have outputted countLeft * grpCountRight
                //tuples. Check if they turned out to be duplicates
                size_t maxSize = countLeft * grpCountRight;
                size_t diff = output->getNRows() - addedSoFar;
                if (!filterDuplEnabled && leftActive && diff < maxSize / 10) {
                    countDuplicatedJoins++;
                    if (retainUnique &&
                            countDuplicatedJoins >= N_ATTEMPTS_ENABLE_DUPL_DEL) {
                        //It seems that the join is producing a huge number
                        //of duplicates. Activate a pre-filtering to prevent
                        //the output of many duplicated derivations
                        if (copyVarPosLeft.size() == 1 &&
                                copyVarPosRight.size() == 1) {
                            filterDuplEnabled = true;
                            LOG(DEBUGL) << "Enabling advanced duplicate "
                                "removal technique";
                            //Populate the two sides of the joins that produce
                            //duplicate derivations
                            duplManager = std::unique_ptr<DuplManager>(
                                    new DuplManager(output.get()));
                        }
                    }
                }
                addedSoFar = output->getNRows();

                countLeft = -1;
                leftGroupDuplicate = false;
                //LeftItr is already pointing to the next element ...
                if (!leftActive) {
                    break;
                }
                res = TGSegmentItr::cmp(itrLeft.get(), itrRight.get(),
                        joinVarsPos);
            } else {
                grpCountRight++;
            }
        }
    }
#if DEBUG
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
                if (shouldTrackProvenance()) {
                    currentrow[sizerow] = itrLeft->getNodeId();
                    if (!copyOnlyLeftNode)
                        currentrow[sizerow + 1] = 0;
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
        if (shouldTrackProvenance()) {
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
        inputLeft = inputLeft->sortBy(fields1);
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
        auto segRight = processAtom_EDB(l, positions);
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
                    if (shouldTrackProvenance()) {
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
