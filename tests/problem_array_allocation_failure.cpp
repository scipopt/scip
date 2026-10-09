// Fail problem-array growth, then release the rejected object and the partial model.
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>
#include <dlfcn.h>
#include <cassert>
#include <cstring>
#include <string>
#include <sys/resource.h>

static bool fail_growth = false;
static bool injected = false;

// Interpose only this executable's calls, including libscip's; no production fault injection.
extern "C" void* BMSreallocMemoryArray_call(void* ptr, size_t count, size_t size,
                                          const char* file, int line) {
    using Function = void* (*)(void*, size_t, size_t, const char*, int);
    static auto original = reinterpret_cast<Function>(dlsym(RTLD_NEXT, "BMSreallocMemoryArray_call"));
    assert(original);
    if (fail_growth && std::strstr(file, "/scip/prob.c")) {
        fail_growth = false;
        injected = true;
        return nullptr;
    }
    return original(ptr, count, size, file, line);
}

static void check(bool constraints) {
    SCIP* scip = nullptr;
    assert(SCIPcreate(&scip) == SCIP_OKAY);
    assert(SCIPincludeDefaultPlugins(scip) == SCIP_OKAY);
    assert(SCIPcreateProbBasic(scip, "allocation_failure") == SCIP_OKAY);
    SCIP_VAR* x = nullptr;
    if (constraints) {
        assert(SCIPcreateVarBasic(scip, &x, "x", 0, 1, 0, SCIP_VARTYPE_CONTINUOUS) == SCIP_OKAY);
        assert(SCIPaddVar(scip, x) == SCIP_OKAY);
    }
    injected = false;
    for (int i = 0; i < 256 && !injected; ++i) {
        fail_growth = i >= 8; // Ensure the array already owns entries when resizing fails.
        const auto name = std::to_string(i);
        SCIP_RETCODE code;
        if (constraints) {
            SCIP_CONS* cons = nullptr;
            double coefficient = 1;
            assert(SCIPcreateConsBasicLinear(scip, &cons, name.c_str(), 1, &x, &coefficient, 0, 1) == SCIP_OKAY);
            code = SCIPaddCons(scip, cons);
            if (injected) {
                assert(SCIPconsGetNUses(cons) == 1);
                assert(!SCIPconsIsAdded(cons));
            }
            assert(SCIPreleaseCons(scip, &cons) == SCIP_OKAY);
        } else {
            SCIP_VAR* var = nullptr;
            assert(SCIPcreateVarBasic(scip, &var, name.c_str(), 0, 1, 0, SCIP_VARTYPE_CONTINUOUS) == SCIP_OKAY);
            code = SCIPaddVar(scip, var);
            if (injected) assert(SCIPvarGetNUses(var) == 1);
            assert(SCIPreleaseVar(scip, &var) == SCIP_OKAY);
        }
        assert(code == (injected ? SCIP_NOMEMORY : SCIP_OKAY));
    }
    assert(injected);
    fail_growth = false;
    if (x) assert(SCIPreleaseVar(scip, &x) == SCIP_OKAY);
    assert(SCIPfree(&scip) == SCIP_OKAY);
    assert(!scip);
}

int main(int argc, char** argv) {
    const rlimit core{0, 0};
    assert(setrlimit(RLIMIT_CORE, &core) == 0);
    // Optional selector permits demonstrating each crash against the unpatched library.
    if (argc == 1 || std::string(argv[1]) == "variables") check(false);
    if (argc == 1 || std::string(argv[1]) == "constraints") check(true);
}
