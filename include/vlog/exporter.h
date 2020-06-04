#include <vlog/chase.h>
#include <vlog/seminaiver.h>
#include <vlog/consts.h>

#include <trident/tree/root.h>

struct _EDBPredicates {
    PredId_t id;
    size_t ruleid;
    int64_t triple[3];
    uint8_t nPosToCopy;
    uint8_t posToCopy[3];
};

class Exporter {
    private:
        std::shared_ptr<Chase> sn;

        void extractTriples(std::vector <uint64_t> &all_s,
                std::vector <uint64_t> &all_p,
                std::vector <uint64_t> &all_o);

        void copyTable(std::vector<uint64_t> &all_s,
                std::vector<uint64_t> &all_p,
                std::vector<uint64_t> &all_o,
                std::vector<_EDBPredicates>::iterator it,
                std::shared_ptr<const FCInternalTable> intTable,
                const int64_t nrows,
                int64_t triple[3]);

    public:
        Exporter(std::shared_ptr<Chase> sn) : sn(sn) {}

        VLIBEXP void generateTridentDiffIndex(std::string outputdir);

        VLIBEXP void generateNTTriples(std::string outputdir, bool decompress);

        VLIBEXP void storeOnFile(std::string path, const PredId_t pred, const bool decompress,
                const int minLevel, const bool csv);

        VLIBEXP void storeOnFiles(std::string path, const bool decompress,
                const int minLevel, const bool csv);
};
