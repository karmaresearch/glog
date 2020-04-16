#include <vlog/text/stringtable_unary.h>

#include <vlog/edb.h>

StringTableUnary::StringTableUnary(PredId_t predid,
        EDBLayer *layer,
        std::string fname,
        std::string param) : StringTable(predid, layer, fname), param(param) {
    buffer1 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
}

bool StringTableUnary::execFunction(const uint64_t t1) {
    bool found1 = layer->getDictText(t1, buffer1.get());
    if (!found1) {
        return false;
    }

    if (fname == "endsWith") {
        auto s = std::string(buffer1.get());
        if (param.size() > s.size()) return false;
        return std::equal(param.rbegin(), param.rend(), s.rbegin());
    }

    if (fname == "isLiteral") {
        auto s = std::string(buffer1.get());
        return s.size() > 1 && s[0] == '\"' && s[s.size() - 1] == '\"';
    }

    LOG(ERRORL) << "(StringTableUnary) Function " << fname << " is unknown";
    throw 10;
}
