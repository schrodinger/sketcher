#include <exception>
#include <string>
#include <vector>

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include "schrodinger/rdkit_extensions/monomer_database.h"

static emscripten::val
string_vector_to_js_array(const std::vector<std::string>& values)
{
    auto array = emscripten::val::array();
    for (unsigned i = 0; i < values.size(); ++i) {
        array.set(i, values[i]);
    }
    return array;
}

static emscripten::val prepare_custom_monomers(const std::string& json)
{
    auto object = emscripten::val::object();
    try {
        auto& db = schrodinger::rdkit_extensions::MonomerDatabase::instance();
        auto result = db.loadMonomersFromJson(json);
        auto data = db.serializeCustomMonomers();
        // Copy out of WASM memory before data is destroyed or transferred.
        object.set("database",
                   emscripten::val::global("Uint8Array")
                       .new_(emscripten::val(emscripten::typed_memory_view(
                           data.size(), data.data()))));
        object.set("succeeded", string_vector_to_js_array(result.first));
        object.set("failed", string_vector_to_js_array(result.second));
    } catch (const std::exception& e) {
        object.set("error", std::string(e.what()));
    }
    return object;
}

EMSCRIPTEN_BINDINGS(monomer_loader)
{
    emscripten::function("prepare_custom_monomers", &prepare_custom_monomers);
}
