// BaumPlane_MEF: ponte entre o núcleo C++ (baumplane_core.hpp) e o navegador (WebAssembly).
// Compilar (sem Emscripten, só clang + wasm-ld):
//   clang++ --target=wasm32 -O2 -std=c++17 -nostdlib -fno-exceptions -fno-rtti -fno-threadsafe-statics \
//           -Wl,--no-entry -Wl,--allow-undefined -Wl,--export-dynamic baumplane_wasm.cpp -o baumplane_mef.wasm
// pow, cos, sin e cbrt vêm do Math do JavaScript (importados de "env"); sqrt e fabs são instruções do wasm.
using bp_real = double;
extern "C" double js_pow(double, double);
extern "C" double js_cos(double);
extern "C" double js_sin(double);
extern "C" double js_cbrt(double);
inline double bp_pow(double a, double b) { return js_pow(a, b); }
inline double bp_cos(double a) { return js_cos(a); }
inline double bp_sin(double a) { return js_sin(a); }
inline double bp_cbrt(double a) { return js_cbrt(a); }
inline double bp_sqrt(double a) { return __builtin_sqrt(a); }
inline double bp_fabs(double a) { return __builtin_fabs(a); }
// sem libc: o compilador gera chamadas a memcpy/memset ao copiar structs
extern "C" void* memcpy(void* d, const void* s, unsigned long n) { auto a = (char*)d; auto b = (const char*)s; while (n--) *a++ = *b++; return d; }
extern "C" void* memset(void* d, int v, unsigned long n) { auto a = (char*)d; while (n--) *a++ = (char)v; return d; }
#include "baumplane_core.hpp"

static double P[bp::NPAR];
static double OUT[16 + 5 * 8 + 1 + 2 * 64];
static bp::Caso C;

#define EXPORT extern "C" __attribute__((visibility("default")))
EXPORT double* bp_params() { for (int i = 0; i < bp::NPAR; i++) if (P[i] == 0) P[i] = bp::PADRAO[i]; return P; }
EXPORT double* bp_out() { return OUT; }
EXPORT double bp_T(double x) { return C.s.T(x); }
// modo 0: estagnação; 1: operação com a placa a Tref. Saída em OUT (ver baumplane_mef.js).
EXPORT int bp_run(int modo, double Tref, int nmax) {
    bool ok = true;
    C = bp::montar(P, 0);
    if (modo == 1) ok = bp::operacao(P, Tref, C);
    const bp::Analitica& s = C.s;
    double v[16] = {C.Tp, s.Tmax(), C.hg, s.k1, s.he, -s.q1(), s.q3(), C.qExtr, C.S, s.L1, s.L2, s.L3, ok ? 1.0 : 0.0, 0, 0, 0};
    bp::Linha L[8];
    int nl = bp::convergencia(s, L, nmax);
    v[13] = nl;
    for (int i = 0; i < 16; i++) OUT[i] = v[i];
    for (int k = 0; k < nl; k++) {
        double* r = OUT + 16 + 5 * k;
        r[0] = L[k].n; r[1] = L[k].N; r[2] = L[k].h; r[3] = L[k].nodal; r[4] = L[k].Er;
    }
    // nós do MEF com n = 1, para o gráfico do perfil
    static double x[bp::NODES], T[bp::NODES];
    int N = bp::malha(s.L1, s.L2, s.L3, 1, x);
    bp::mef(s, x, N, T);
    double* q = OUT + 56;
    q[0] = N;
    for (int i = 0; i < N; i++) { q[1 + 2 * i] = x[i]; q[2 + 2 * i] = T[i]; }
    return nl;
}
