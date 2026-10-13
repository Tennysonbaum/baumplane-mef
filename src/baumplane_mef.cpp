// BaumPlane_colector: programa de linha de comando do núcleo MEF (baumplane_core.hpp).
// Mesma física do app e do script fcc220_reprodutibilidade.py (Apêndices A e C do TCC).
//
// Compilar:  g++ -O2 -std=c++17 baumplane_mef.cpp -o baumplane_mef
// Precisão estendida (long double, estende a tabela de convergência como no TCC):
//            g++ -O2 -std=c++17 -DESTENDIDA baumplane_mef.cpp -o baumplane_mef_ld
// Rodar:     ./baumplane_mef                    (FCC220-2V, estagnação)
//            ./baumplane_mef caso=op            (operação com a placa a 60 °C)
//            ./baumplane_mef eiso=40 gap=16 eg=4   (troca qualquer parâmetro: chave=valor)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifdef ESTENDIDA
using bp_real = long double;
#define RF "Lf"
#define RE "Le"
#else
using bp_real = double;
#define RF "f"
#define RE "e"
#endif
inline bp_real bp_pow(bp_real a, bp_real b) { return std::pow(a, b); }
inline bp_real bp_cos(bp_real a) { return std::cos(a); }
inline bp_real bp_sin(bp_real a) { return std::sin(a); }
inline bp_real bp_cbrt(bp_real a) { return std::cbrt(a); }
inline bp_real bp_sqrt(bp_real a) { return std::sqrt(a); }
inline bp_real bp_fabs(bp_real a) { return std::fabs(a); }
#include "baumplane_core.hpp"
using bp::real;

int main(int argc, char** argv) {
    real P[bp::NPAR];
    for (int i = 0; i < bp::NPAR; i++) P[i] = bp::PADRAO[i];
    std::string caso = "stag";
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i]; size_t eq = a.find('=');
        if (eq == std::string::npos) { std::fprintf(stderr, "argumento invalido: %s (use chave=valor)\n", argv[i]); return 1; }
        std::string k = a.substr(0, eq), v = a.substr(eq + 1);
        if (k == "caso") { caso = v; continue; }
        int j = 0;
        while (j < bp::NPAR && k != bp::NOMES[j]) j++;
        if (j == bp::NPAR) { std::fprintf(stderr, "parametro desconhecido: %s\n", k.c_str()); return 1; }
        P[j] = std::strtold(v.c_str(), nullptr);
    }
    bp::Caso C = bp::montar(P, 0);
    if (caso == "op" && !bp::operacao(P, 60, C)) std::printf("Em estagnacao a placa nao chega a 60 C; mostrando estagnacao.\n");
    const bp::Analitica& s = C.s;
    real qf = -s.q1(), qb = s.q3(), dT = C.Tp - P[bp::TA];
    std::printf("BaumPlane_colector: MEF 1D (%s, %s)\n", caso == "op" ? "operacao, placa a 60 C" : "estagnacao",
                sizeof(real) > 8 ? "long double" : "double");
    std::printf("  placa T(L1)         = %.3" RF " C   (Tmax no absorvedor %.6" RF " C)\n", C.Tp, s.Tmax());
    std::printf("  h camara            = %.3" RF " W/m2K, k1 equivalente = %.4" RF " W/mK, h externo = %.2" RF " W/m2K\n", C.hg, s.k1, s.he);
    std::printf("  balanco             : frente %.2" RF " + fundo %.2" RF " + fluido %.2" RF " = S %.2" RF " W/m2 (residuo %.1" RE ")\n",
                qf, qb, C.qExtr, C.S, std::fabs(qf + qb + C.qExtr - C.S));
    std::printf("  U frente + fundo    = %.4" RF " + %.4" RF " = %.4" RF " W/m2K\n\n", qf / dT, qb / dT, (qf + qb) / dT);

    bp::Linha L[8];
    int nl = bp::convergencia(s, L);
    std::printf("   n    nos   h_II (mm)   erro nodal max (K)   E_r (L2)      p\n");
    bool piso = false;                       // depois que o arredondamento aparece, não volta
    for (int k = 0; k < nl; k++) {
        std::printf("  %3d  %5d   %9.5" RF "   %18.2" RE "   %9.2" RE, L[k].n, L[k].N, L[k].h, L[k].nodal, L[k].Er);
        if (k) { real p = std::log2(L[k - 1].Er / L[k].Er); piso = piso || std::fabs(p - 2) >= 0.15; std::printf("   %5.2" RF "%s", p, piso ? "  (arredondamento)" : ""); }
        else std::printf("       -");
        std::printf("\n");
    }
    return 0;
}
