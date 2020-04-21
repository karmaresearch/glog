#include <vlog/text/stringtable_unary.h>

#include <vlog/edb.h>

StringTableUnary::StringTableUnary(PredId_t predid,
        EDBLayer *layer,
        std::string fname,
        std::string param) : StringTable(predid, layer, fname), param(param) {
    buffer1 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
    if (fname == "maxLen") {
        param1_int = atoi(param.c_str());
    }
    if (fname == "containsNoChar") {
        param1_char = param.c_str()[0];
    }
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

    if (fname == "maxLen") {
        auto len = strlen(buffer1.get());
        bool outcome = len <= param1_int;
        return outcome;
    }

    if (fname == "containsNoChar") {
        auto pos = strchr (buffer1.get(), param1_char);
        bool outcome = pos == NULL;
        return outcome;
    }

    LOG(ERRORL) << "(StringTableUnary) Function " << fname << " is unknown";
    throw 10;
}
