#ifndef _PROB_GBCHASE_H
#define _PROB_GBCHASE_H

#include <glog/gbchase.h>

#include <vlog/edb.h>
#include <vlog/concepts.h>

class ProbGBChase : public GBChase
{
    protected:
        bool executeRule(GBRuleInput &node, bool cleanDuplicates = true);

    public:
        ProbGBChase(EDBLayer &layer, Program &program);
};
#endif
