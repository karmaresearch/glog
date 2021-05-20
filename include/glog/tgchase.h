#ifndef _TG_CHASE_H
#define _TG_CHASE_H

#include <vlog/gbchase/gbgraph.h>
#include <vlog/gbchase/gbchase.h>

#include <chrono>

typedef std::vector<std::string> TGChase_NodeIDs;

class TGChase_RuleIO {
    public:
        TGChase_NodeIDs outNodes;
        std::vector<size_t> outSizes;
        std::vector<TGChase_NodeIDs> inNodes;

        TGChase_RuleIO(const std::string &outNode,
                const std::vector<std::string> &inNodes,
                bool removeEDBNodes = true) {
            outNodes.push_back(outNode);
            for(auto &n : inNodes) {
                if (removeEDBNodes && n.rfind("EDB", 0) == 0) {
                    continue;
                }
                std::vector<std::string> v;
                v.push_back(n);
                this->inNodes.push_back(v);
            }
        }

};

class TGChase : public GBChase {
    private:
        std::unordered_map<std::string, size_t> mapExternalInternalIDs;

    public:
        VLIBEXP TGChase(EDBLayer &layer, Program *program);

        virtual void run() = 0;

        void executeRule(const size_t ruleID, TGChase_RuleIO &io);
};

class TGChaseStatic : public TGChase {
    private:
        std::string tgfile;

    public:
        VLIBEXP TGChaseStatic(EDBLayer &layer, Program *program,
                std::string tgfile) : TGChase(layer, program), tgfile(tgfile) {
        }

        VLIBEXP virtual void run();
};

#endif
