#ifndef _GB_SEGMENT_INSERTER_H
#define _GB_SEGMENT_INSERTER_H

#include <vlog/gbchase/gbsegment.h>

#include <vlog/term.h>
#include <vlog/column.h>
#include <vlog/segment.h>

#include <cstddef>
#include <memory>

class GBSegmentInserter {
    public:
        static std::unique_ptr<GBSegmentInserter> getInserter(size_t card);

        virtual bool isEmpty() const = 0;

        virtual size_t getNRows() const = 0;

        virtual void addRow(Term_t *row) = 0;

        virtual std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance) = 0;

        virtual void postprocessJoin(std::vector<std::shared_ptr<Column>>
                &intermediateResultsNodes) = 0;

        virtual ~GBSegmentInserter() {}
};

template<class K>
class GBSegmentInserterImpl : public GBSegmentInserter {
    protected:
        const size_t card;
        K tuples;

    public:
        GBSegmentInserterImpl(size_t card) : card(card) {
        }

        bool isEmpty() const {
            return tuples.empty();
        }

        virtual size_t getNRows() const {
            return tuples.size();
        }

        void postprocessJoin(std::vector<std::shared_ptr<Column>>
                &intermediateResultsNodes) {
            LOG(ERRORL) << "Should not occur";
            throw 10;
        }
};

class GBSegmentInserterUnary : public GBSegmentInserterImpl<std::vector<Term_t>>
{
    public:
        GBSegmentInserterUnary() : GBSegmentInserterImpl(1) { }

        void addRow(Term_t *row) {
            tuples.push_back(row[0]);
        }

        std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance) {
            assert(trackProvenance == false);
            return std::shared_ptr<const TGSegment>(
                    new UnaryTGSegment(tuples, nodeId, isSorted, sortedField));
        }
};

class GBSegmentInserterBinary : public GBSegmentInserterImpl<
                                std::vector<std::pair<Term_t, Term_t>>>
{
    private:
        bool secondFieldConstant;

    public:
        GBSegmentInserterBinary() : GBSegmentInserterImpl(1),
        secondFieldConstant(true) { }

        void addRow(Term_t *row) {
            if (secondFieldConstant && !tuples.empty() &&
                    tuples.back().second != row[1]) {
                secondFieldConstant = false;
            }
            tuples.push_back(std::make_pair(row[0], row[1]));
        }

        std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance) {
            if (trackProvenance) {
                if (secondFieldConstant) {
                    std::vector<Term_t> t;
                    for(auto tuple : tuples) {
                        t.push_back(tuple.first);
                    }
                    return std::shared_ptr<const TGSegment>(
                            new UnaryWithConstProvTGSegment(
                                t, nodeId, isSorted, sortedField));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new UnaryWithProvTGSegment(tuples, ~0ul, isSorted,
                                sortedField));
                }
            } else {
                return std::shared_ptr<const TGSegment>(
                        new BinaryTGSegment(tuples, nodeId, isSorted, sortedField));
            }
        }
};

class GBSegmentInserterNAry : public GBSegmentInserter
{
    private:
        SegmentInserter ins;
        const size_t card;
        std::shared_ptr<const Segment> seg;
        bool inserterIsClosed;

    public:
        GBSegmentInserterNAry(size_t card) : ins(card), card(card),
        inserterIsClosed(false) { }

        bool isEmpty() const  {
            return ins.isEmpty();
        }

        size_t getNRows() const {
            return ins.getNRows();
        }

        void addRow(Term_t *row) {
            assert(!inserterIsClosed);
            ins.addRow(row);
        }

        std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance) {
            if (!inserterIsClosed) {
                seg = ins.getSegment();
                inserterIsClosed = true;
            }
            if (card == 3 && trackProvenance) {
                auto &col1 = seg->getColumn(0)->getVectorRef();
                auto &col2 = seg->getColumn(1)->getVectorRef();
                bool constantNodeVal = seg->getColumn(2)->isConstant();
                size_t nrows = col1.size();
                if (constantNodeVal) {
                    std::vector<std::pair<Term_t, Term_t>> out(nrows);
                    for(size_t i = 0; i < nrows; ++i) {
                        out[i].first = col1[i];
                        out[i].second = col2[i];
                    }
                    return std::shared_ptr<const TGSegment>(
                            new BinaryWithConstProvTGSegment(
                                out, seg->getColumn(2)->first(), isSorted, sortedField));
                } else {
                    std::vector<BinWithProv> out(nrows);
                    auto itrNode = seg->getColumn(2)->getReader();
                    for(size_t i = 0; i < nrows; ++i) {
                        out[i].first = col1[i];
                        out[i].second = col2[i];
                        if (!itrNode->hasNext()) {
                            throw 10;
                        }
                        out[i].node = itrNode->next();
                    }
                    return std::shared_ptr<const TGSegment>(
                            new BinaryWithProvTGSegment(
                                out, ~0ul, isSorted, sortedField));
                }
            } else {
                std::vector<std::shared_ptr<Column>> columns;
                for(int i = 0; i < card; ++i) {
                    columns.push_back(seg->getColumn(i));
                }
                if (trackProvenance) {
                    return std::shared_ptr<const TGSegment>(
                            new TGSegmentLegacy(columns, seg->getNRows(), isSorted,
                                sortedField, true));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new TGSegmentLegacy(columns, seg->getNRows(), isSorted,
                                sortedField, false));
                }
            }
        }

        void postprocessJoin(std::vector<std::shared_ptr<Column>>
                &intermediateResultsNodes) {
            assert(!inserterIsClosed);
            inserterIsClosed = true;
            seg = ins.getSegment();

            //Take out the last two columns
            auto ncolumns = seg->getNColumns();
            std::vector<std::shared_ptr<Column>> columns;
            for(int i = 0; i < ncolumns - 2; ++i) {
                columns.push_back(seg->getColumn(i));
            }
            //For now, always add one extra column
            CompressedColumnBlock b(0, 1, columns[0]->size());
            std::vector<CompressedColumnBlock> blocks;
            blocks.push_back(b);
            columns.push_back(std::shared_ptr<Column>(
                        new CompressedColumn(blocks, columns[0]->size())));

            //Save the columns with the mappings to the nodes
            auto col1 = seg->getColumn(ncolumns - 2);
            auto col2 = seg->getColumn(ncolumns - 1);
            intermediateResultsNodes.push_back(col1);
            intermediateResultsNodes.push_back(col2);

            //Create new intermediate results
            seg = std::shared_ptr<const Segment>(new Segment(columns.size(), columns));

        }

};

#endif
