#include <vlog/text/elastictable.h>
#include <vlog/text/elasticitr.h>

#include <kognac/logs.h>
#include <trident/utils/httpclient.h>
#include <trident/utils/json.h>

ElasticTable::ElasticTable(PredId_t predid,
        EDBLayer *layer,
        std::string baserel,
        std::string basehost,
        std::string baseport,
        std::string basepath,
        std::string startRange) : predid(predid), layer(layer), baserel(baserel),
    basehost(basehost), baseport(stoi(baseport)),
    basepath(basepath), startRange(stol(startRange)) {
        //get nterms from the baserel
        auto dictPredId = layer->getPredID(baserel);
        dictTable = layer->getEDBTable(dictPredId);
        nterms = dictTable->getNTerms();
    }

void ElasticTable::query(QSQQuery *query, TupleTable *outputTable,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t ElasticTable::getCardinality(const Literal &query) {
    return 1;
}

size_t ElasticTable::getCardinalityColumn(const Literal &query,
        uint8_t posColumn) {
    return 1;
}

bool ElasticTable::isQueryAllowed(const Literal &query) {
    auto t1 = query.getTermAtPos(0);
    return !t1.isVariable();
}

size_t ElasticTable::estimateCardinality(const Literal &query) {
    return getCardinality(query);
}

bool ElasticTable::isEmpty(const Literal &query,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *ElasticTable::getIterator(const Literal &query) {
    auto t1 = query.getTermAtPos(0);
    assert(!t1.isVariable());
    return new ElasticIterator(predid, t1.getValue(), t1.getValue() + startRange);
}

EDBIterator *ElasticTable::getSortedIterator(const Literal &query,
        const std::vector<uint8_t> &fields) {
    return getIterator(query);
}

bool ElasticTable::getDictNumber(const char *text, const size_t sizeText,
        uint64_t &id) {
    return false;
}

bool ElasticTable::getDictText(const uint64_t id, char *text) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool ElasticTable::getDictText(const uint64_t id, std::string &text) {
    auto origid = id;
    if (id >= startRange && id < startRange + nterms) {
        origid -= startRange;
    }
    auto resp = dictTable->getDictText(origid, text);
    if (resp && id >= startRange && id < startRange + nterms) {
        //Prepare the elasticsearch request
        std::string params = "{ \"query\": { \"match\": { \"_id\": \"" + text + "\"}}}";
        std::string headers = "";
        std::string contenttype = "application/json";
        std::string response;
        HttpClient client(basehost, baseport);
        client.connect();
        resp = client.post(basepath, params, headers, response, contenttype);
        bool ok = false;
        if (resp) {
            std::cout << response << std::endl;
            JSON r;
            JSON::read(response, r);
            if (r.containsChild("hits")) {
                auto h = r.getChild("hits");
		if (h.containsChild("hits")) {
		auto h2 = h.getChild("hits");
                auto &hits = h2.getListChildren();
                if (hits.size() > 0) {
                    JSON hit = hits[0];
                    if (hit.containsChild("_source")) {
			    auto fields = hit.getChild("_source");
			    if (fields.contains("rdfs_label")) {
                        text = hit.get("rdfs_label");
                        ok = true;
                    }
		}
                }
		}
            }
        }
        if (!ok) {
            text = "LABEL_NOT_FOUND";
        }
    }
    return resp;
}

uint64_t ElasticTable::getNTerms() {
    return nterms;
}

uint64_t ElasticTable::getSize() {
    return nterms;
}

void ElasticTable::releaseIterator(EDBIterator *itr) {
    delete itr;
}

ElasticTable::~ElasticTable() {
}
