#include <vlog/text/stringtable_binary.h>

#include <vlog/edb.h>

#include <cstring>

StringTableBinary::StringTableBinary(PredId_t predid,
        EDBLayer *layer,
        std::string fname) : StringTable(predid, layer, fname) {
    buffer1 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
    buffer2 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
}

bool StringTableBinary::execFunction(const uint64_t t1, const uint64_t t2) {
    bool found1 = layer->getDictText(t1, buffer1.get());
    bool found2 = layer->getDictText(t2, buffer2.get());
    if (!found1 || !found2) {
        return false;
    }

    if (fname == "containedIn") {
        auto resp = std::string(buffer2.get()).find(std::string(buffer1.get()));
        return resp != std::string::npos;
    } else if (fname == "equal") {
        return strcmp(buffer1.get(), buffer2.get()) == 0;
    }

    LOG(ERRORL) << "(StringTableBinary) Function " << fname << " is unknown";
    throw 10;
}
