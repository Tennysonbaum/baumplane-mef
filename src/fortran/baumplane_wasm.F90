! BaumPlane_MEF: o núcleo Fortran (módulo baumplane de baumplane_mef.F90) exportado para WebAssembly.
! O mesmo código-fonte do programa de linha de comando; aqui só se acrescentam funções bind(c) para o JavaScript.
! Compilar (flang-new 18 + clang 18; o flang 18 ainda não gera wasm32, então gera-se o LLVM IR para i386
! e troca-se o alvo para wasm32, o que vale porque o módulo não chama a biblioteca de execução do Fortran):
!   flang-new -fc1 -triple i386-unknown-linux-gnu -emit-llvm -O2 -cpp baumplane_wasm.F90 -o bw.ll
!   sed -e "s/^target datalayout = .*/target datalayout = \"e-m:e-p:32:32-p10:8:8-p20:8:8-i64:64-n32:64-S128-ni:1:10:20\"/" \
!       -e "s/^target triple = .*/target triple = \"wasm32-unknown-unknown\"/" bw.ll > bw32.ll
!   clang --target=wasm32 -O2 -nostdlib -fno-builtin -Wl,--no-entry -Wl,--allow-undefined -Wl,-z,stack-size=1048576 \
!         -Wl,--export=bp_set -Wl,--export=bp_get -Wl,--export=bp_T -Wl,--export=bp_run bw32.ll shim.c \
!         libclang_rt.builtins-wasm32.a -o baumplane_mef_f.wasm
! O JavaScript fornece pow, cos, sin e round. Para a precisão quádrupla (baumplane_mef_fq.wasm), acrescentar
! -DESTENDIDA -DMATEMATICA_PROPRIA na linha do flang: o real(16) vira aritmética em software (as rotinas
! __addtf3, __multf3... de libclang_rt.builtins-wasm32.a, pacote libclang-rt-18-dev-wasm32) e potência, seno,
! cosseno e raiz saem das séries do próprio módulo.
#define SO_MODULO
#include "baumplane_mef.F90"

module baumplane_js
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use baumplane
  implicit none
  real(wp), save :: P(NPAR) = PADRAO
  real(wp), save :: OUT(200) = 0
  type(caso), save :: C
  real(wp), save :: gp(3), gw(3)                      ! pontos e pesos de Gauss (3 pontos)
contains
  subroutine bp_set(i, v) bind(c, name='bp_set')            ! parâmetro i (1..16) = v
    integer(c_int), value :: i
    real(c_double), value :: v
    real(c_double) :: e, m
    integer :: k
    if (i < 1 .or. i > NPAR) return
    P(i) = v
    e = 1
    do k = 0, 9                                ! decimal curto (0.4, 237, 0.88): refaz a divisão na precisão do
      m = anint(v*e)                           ! núcleo, como faria a leitura de "0.4" pela linha de comando
      if (abs(v*e - m) <= 1e-9_c_double*max(1._c_double, abs(m))) then
        P(i) = real(m, wp)/real(e, wp); return
      end if
      e = e*10
    end do
  end subroutine
  real(c_double) function bp_get(i) bind(c, name='bp_get')  ! saída i (1..200), ver baumplane_mef.js
    integer(c_int), value :: i
    bp_get = 0
    if (i >= 1 .and. i <= 200) bp_get = OUT(i)
  end function
  real(c_double) function bp_T(xv) bind(c, name='bp_T')     ! solução analítica em x (m)
    real(c_double), value :: xv
    bp_T = T_an(C%s, real(xv, wp))
  end function
  integer(c_int) function bp_run(modo, Tref, nmax) bind(c, name='bp_run')
    integer(c_int), value :: modo, nmax
    real(c_double), value :: Tref
    real(wp) :: lo, hi_, m, x(NODES), T(NODES), nrm, nod, e2
    integer :: i, e, n, N_, k
    logical :: ok
    ok = .true.
    call montar(P, 0._wp, C)
    if (modo == 1) then
      if (C%Tp <= Tref) then
        ok = .false.
      else
        lo = 0; hi_ = C%Sabs
        do i = 1, 80
          m = (lo + hi_)/2; call montar(P, m, C)
          if (C%Tp > Tref) then
            lo = m
          else
            hi_ = m
          end if
        end do
      end if
    end if
    associate (s => C%s)
      OUT(1) = C%Tp; OUT(2) = s%k2*s%c(3)*s%c(3)/(2*s%q) + s%c(4); OUT(3) = C%hg; OUT(4) = s%k1; OUT(5) = s%he
      OUT(6) = s%k1*s%c(1); OUT(7) = -s%k3*s%c(5); OUT(8) = C%qExtr; OUT(9) = C%Sabs
      OUT(10) = s%L1; OUT(11) = s%L2; OUT(12) = s%L3; OUT(13) = merge(1._wp, 0._wp, ok)
      gp(1) = -raiz(0.6_wp); gp(2) = 0; gp(3) = raiz(0.6_wp); gw(1) = 5._wp/9; gw(2) = 8._wp/9; gw(3) = 5._wp/9
      nrm = raiz(q2(0._wp, s%L1) + q2(s%L1, s%L2) + q2(s%L2, s%L3))
      n = 1; k = 0
      do while (n <= nmax .and. n <= NMAX)
        N_ = malha(s%L1, s%L2, s%L3, n, x)
        call mef(s, x, N_, T)
        nod = 0; e2 = 0
        do i = 1, N_
          nod = max(nod, abs(T(i) - T_an(s, x(i))))
        end do
        do e = 1, N_ - 1
          e2 = e2 + qe(x(e), x(e + 1), T(e), T(e + 1))
        end do
        OUT(17 + 5*k) = n; OUT(18 + 5*k) = N_; OUT(19 + 5*k) = (s%L2 - s%L1)/n*1000
        OUT(20 + 5*k) = nod; OUT(21 + 5*k) = raiz(e2)/nrm
        k = k + 1; n = 2*n
      end do
      OUT(14) = k
      N_ = malha(s%L1, s%L2, s%L3, 1, x)                     ! nós com n = 1, para o gráfico
      call mef(s, x, N_, T)
      OUT(57) = N_
      do i = 1, N_
        OUT(56 + 2*i) = x(i); OUT(57 + 2*i) = T(i)
      end do
    end associate
    bp_run = k
  end function
  real(wp) function q2(a, b)
    real(wp), intent(in) :: a, b
    integer :: j
    real(wp) :: z, t
    q2 = 0
    do j = 1, 3
      z = (b - a)/2*gp(j) + (a + b)/2; t = T_an(C%s, z); q2 = q2 + gw(j)*t*t
    end do
    q2 = (b - a)/2*q2
  end function
  real(wp) function qe(a, b, ta_, tb)
    real(wp), intent(in) :: a, b, ta_, tb
    integer :: j
    real(wp) :: z, d
    qe = 0
    do j = 1, 3
      z = (b - a)/2*gp(j) + (a + b)/2; d = ta_ + (tb - ta_)*(z - a)/(b - a) - T_an(C%s, z); qe = qe + gw(j)*d*d
    end do
    qe = (b - a)/2*qe
  end function
end module
