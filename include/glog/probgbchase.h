#ifndef _PROB_GBCHASE_H
#define _PROB_GBCHASE_H

#include <glog/gbchase.h>

#include <vlog/edb.h>
#include <vlog/concepts.h>

class ProbGBChase : public GBChase
{
    private:
        const bool optDelProofsStaticAnalysis;

    protected:
        void prepareRuleExecutionPlans(
                const size_t &ruleIdx,
                const size_t prevstep,
                const size_t step,
                std::vector<GBRuleInput> &newnodes);

        bool executeRule(GBRuleInput &node, bool cleanDuplicates = true);

    public:
        ProbGBChase(EDBLayer &layer, Program &program,
                bool optDelProofsStaticAnalysis);
};
#endif
