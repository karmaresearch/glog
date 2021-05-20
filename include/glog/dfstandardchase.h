#ifndef _DFS_CHASE_H
#define _DFS_CHASE_H

#include <vlog/gbchase/gbchase.h>

class DFStandardChase : public GBChase {
    public:
        size_t executeRulesInStratum(
                const std::vector<size_t> &ruleIdxs,
                const size_t stratumLevel,
                const size_t stepStratum,
                size_t &step);
};

#endif
