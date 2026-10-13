// BaumPlane_colector: núcleo de cálculo (modelo por regiões + MEF linear 1D), Apêndices A e C do TCC.
// Não usa a biblioteca padrão: o mesmo código compila nativo (g++) e para WebAssembly (clang).
// Quem inclui define `bp_real` e as funções bp_pow, bp_cos, bp_sin, bp_cbrt, bp_sqrt, bp_fabs, bp_log2.
#pragma once
namespace bp {
using real = bp_real;
#define RL(x) ((real)x##L)   // literal na precisão de real, sem promover a long double
const real SIG = RL(5.670374e-8), PI = RL(3.14159265358979323846);
const int NMAX = 64, NODES = 13 * NMAX + 1;          // malha 4n/n/8n até n = 64

// Parâmetros na ordem usada pelo app e pela linha de comando (espessuras em mm)
enum { G, TA, BETA, HW, HI, EG, KG, EPSG, TAU, GAP, EABS, KABS, ALPHA, EPS, EISO, KISO, NPAR };
static const char* NOMES[NPAR] = {"G", "Ta", "beta", "hw", "hi", "eg", "kg", "epsg", "tau", "gap", "eabs", "kabs", "alpha", "eps", "eiso", "kiso"};
static const real PADRAO[NPAR] = {800, 25, 45, 10, 5, RL(3.2), RL(1.0), RL(0.88), RL(0.91), 25, RL(0.4), 237, RL(0.95), RL(0.05), 25, RL(0.040)};  // FCC220-2V

inline real mx(real a, real b) { return a > b ? a : b; }

// ---------- correlações (Duffie e Beckman, 2013, cap. 3) ----------
inline real hGap(real Tp, real Tg, real L, real b) {     // Hollands, eq. 3.11.4
    real TmK = (Tp + Tg) / 2 + RL(273.15);
    real k = RL(0.0241) * bp_pow(TmK / RL(273.15), RL(0.81)), nu = RL(1.33e-5) * bp_pow(TmK / RL(273.15), RL(1.75)), a = nu / RL(0.71);
    real dT = mx(bp_fabs(Tp - Tg), RL(1e-6));
    real Ra = RL(9.81) / TmK * dT * L * L * L / (nu * a), c = bp_cos(b * PI / 180);
    real t1 = mx(0, 1 - 1708 / (Ra * c));
    real Nu = 1 + RL(1.44) * (1 - 1708 * bp_pow(bp_sin(RL(1.8) * b * PI / 180), RL(1.6)) / (Ra * c)) * t1
              + mx(0, bp_cbrt(Ra * c / 5830) - 1);
    return Nu * k / L;
}
inline real hRad(real T1, real T2, real e1, real e2) {   // placas paralelas, eq. 3.10.2
    real a = T1 + RL(273.15), b = T2 + RL(273.15);
    return SIG * (a * a + b * b) * (a + b) / (1 / e1 + 1 / e2 - 1);
}

// ---------- solução analítica por regiões: sistema 6×6 (Apêndice A) ----------
struct Analitica {
    real k1, k2, k3, L1, L2, L3, he, hi, Tinf, q, c[6];   // c = A1 B1 A2 B2 A3 B3
    real T(real x) const {
        if (x <= L1) return c[0] * x + c[1];
        if (x <= L2) return -q / (2 * k2) * x * x + c[2] * x + c[3];
        return c[4] * x + c[5];
    }
    real q1() const { return -k1 * c[0]; }    // < 0: sai pela frente
    real q3() const { return -k3 * c[4]; }    // > 0: sai pelo fundo
    real Tmax() const { return k2 * c[2] * c[2] / (2 * q) + c[3]; }
};
inline void resolver6(real M[6][6], real b[6], real x[6]) {   // Gauss com pivotamento parcial
    for (int c = 0; c < 6; c++) {
        int p = c;
        for (int r = c + 1; r < 6; r++) if (bp_fabs(M[r][c]) > bp_fabs(M[p][c])) p = r;
        for (int k = 0; k < 6; k++) { real t = M[c][k]; M[c][k] = M[p][k]; M[p][k] = t; }
        real t = b[c]; b[c] = b[p]; b[p] = t;
        for (int r = c + 1; r < 6; r++) {
            real f = M[r][c] / M[c][c];
            for (int k = c; k < 6; k++) M[r][k] -= f * M[c][k];
            b[r] -= f * b[c];
        }
    }
    for (int r = 5; r >= 0; r--) {
        real s = b[r];
        for (int k = r + 1; k < 6; k++) s -= M[r][k] * x[k];
        x[r] = s / M[r][r];
    }
}
inline Analitica analitica(real k1, real k2, real k3, real L1, real L2, real L3, real he, real hi, real Tinf, real q) {
    real M[6][6] = {{k1, -he, 0, 0, 0, 0}, {0, 0, 0, 0, -(k3 + hi * L3), -hi}, {L1, 1, -L1, -1, 0, 0},
                    {-k1, 0, k2, 0, 0, 0}, {0, 0, L2, 1, -L2, -1}, {0, 0, -k2, 0, k3, 0}};
    real b[6] = {-he * Tinf, -hi * Tinf, -q * L1 * L1 / (2 * k2), q * L1, q * L2 * L2 / (2 * k2), -q * L2};
    Analitica s = {k1, k2, k3, L1, L2, L3, he, hi, Tinf, q, {0, 0, 0, 0, 0, 0}};
    resolver6(M, b, s.c);
    return s;
}

// ---------- MEF linear (Apêndice C): matriz tridiagonal resolvida por Thomas ----------
inline int malha(real L1, real L2, real L3, int n, real* x) {    // 4n / n / 8n elementos; devolve nº de nós
    int N = 0;
    for (int i = 0; i <= 4 * n; i++) x[N++] = L1 * i / (4 * n);
    for (int i = 1; i <= n; i++) x[N++] = L1 + (L2 - L1) * i / n;
    for (int i = 1; i <= 8 * n; i++) x[N++] = L2 + (L3 - L2) * i / (8 * n);
    return N;
}
inline void mef(const Analitica& s, const real* x, int N, real* T) {
    static real a[NODES], d[NODES], c[NODES], F[NODES];
    for (int i = 0; i < N; i++) a[i] = d[i] = c[i] = F[i] = 0;
    for (int e = 0; e + 1 < N; e++) {
        real xm = (x[e] + x[e + 1]) / 2, Le = x[e + 1] - x[e];
        real k = xm < s.L1 ? s.k1 : xm < s.L2 ? s.k2 : s.k3, qv = (xm > s.L1 && xm < s.L2) ? s.q : 0;
        real st = k / Le, f = qv * Le / 2;
        d[e] += st; d[e + 1] += st; c[e] -= st; a[e + 1] -= st; F[e] += f; F[e + 1] += f;
    }
    d[0] += s.he; F[0] += s.he * s.Tinf;               // Robin na frente
    d[N - 1] += s.hi; F[N - 1] += s.hi * s.Tinf;       // Robin no fundo
    for (int i = 1; i < N; i++) { real m = a[i] / d[i - 1]; d[i] -= m * c[i - 1]; F[i] -= m * F[i - 1]; }
    T[N - 1] = F[N - 1] / d[N - 1];
    for (int i = N - 2; i >= 0; i--) T[i] = (F[i] - c[i] * T[i + 1]) / d[i];
}

// ---------- montagem com ponto fixo dos coeficientes (câmara e radiação para o céu) ----------
struct Caso { Analitica s; real hg, qExtr, S, Tp; };
inline Caso montar(const real* P, real qExtr) {
    real eg = P[EG] / 1000, Lg = P[GAP] / 1000, eabs = P[EABS] / 1000;
    real L1 = eg + Lg, L2 = L1 + eabs, L3 = L2 + P[EISO] / 1000;
    real S = P[G] * RL(1.01) * P[TAU] * P[ALPHA], qv = (S - qExtr) / eabs;
    real Tp = 80, Tg = 40, To = 35, hg = 0;
    Analitica s = {};
    for (int it = 0; it < 200; it++) {
        hg = hGap(Tp, Tg, Lg, P[BETA]) + hRad(Tp, Tg, P[EPS], P[EPSG]);
        real k1 = L1 / (eg / P[KG] + 1 / hg);
        real A = To + RL(273.15), B = P[TA] + RL(273.15);
        real he = P[HW] + P[EPSG] * SIG * (A * A + B * B) * (A + B);   // Robin + radiação linearizada
        s = analitica(k1, P[KABS], P[KISO], L1, L2, L3, he, P[HI], P[TA], qv);
        real nTp = s.T(L1), nTo = s.T(0), nTg = nTo + bp_fabs(s.q1()) * eg / P[KG];
        real dd = mx(mx(bp_fabs(nTp - Tp), bp_fabs(nTg - Tg)), bp_fabs(nTo - To));
        Tp = nTp; Tg = nTg; To = nTo;
        if (dd < RL(1e-12)) break;
    }
    Caso C = {s, hg, qExtr, S, Tp};
    return C;
}
// caso "op": extração que deixa a placa a Tref (bissecção); devolve false se a estagnação não chega lá
inline bool operacao(const real* P, real Tref, Caso& C) {
    C = montar(P, 0);
    if (C.Tp <= Tref) return false;
    real lo = 0, hi = C.S;
    for (int i = 0; i < 80; i++) { real m = (lo + hi) / 2; C = montar(P, m); if (C.Tp > Tref) lo = m; else hi = m; }
    return true;
}

// ---------- convergência: erro nodal máximo e erro relativo L² (Gauss 3 pontos) ----------
struct Linha { int n, N; real h, nodal, Er; };
inline int convergencia(const Analitica& s, Linha* L, int nmax = NMAX) {
    static real x[NODES], T[NODES];
    const real gp[3] = {-bp_sqrt(RL(0.6)), 0, bp_sqrt(RL(0.6))}, gw[3] = {RL(5.0) / 9, RL(8.0) / 9, RL(5.0) / 9};
    auto quad = [&](real a, real b, auto f) { real r = 0; for (int j = 0; j < 3; j++) r += gw[j] * f((b - a) / 2 * gp[j] + (a + b) / 2); return (b - a) / 2 * r; };
    auto T2 = [&](real z) { real t = s.T(z); return t * t; };
    real nrm = bp_sqrt(quad(0, s.L1, T2) + quad(s.L1, s.L2, T2) + quad(s.L2, s.L3, T2));   // exato: T² tem grau ≤ 4
    int k = 0;
    for (int n = 1; n <= nmax; n *= 2, k++) {
        int N = malha(s.L1, s.L2, s.L3, n, x);
        mef(s, x, N, T);
        real nod = 0, e2 = 0;
        for (int i = 0; i < N; i++) nod = mx(nod, bp_fabs(T[i] - s.T(x[i])));
        for (int e = 0; e + 1 < N; e++) {
            real a = x[e], b = x[e + 1], Ta = T[e], Tb = T[e + 1];
            e2 += quad(a, b, [&](real z) { real d = Ta + (Tb - Ta) * (z - a) / (b - a) - s.T(z); return d * d; });
        }
        L[k] = {n, N, (s.L2 - s.L1) / n * 1000, nod, bp_sqrt(e2) / nrm};
    }
    return k;
}
}  // namespace bp
