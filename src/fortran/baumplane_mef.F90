! BaumPlane_colector: MEF 1D em Fortran (modelo por regiões, Apêndices A e C do TCC).
! Mesma física e mesmo algoritmo de baumplane_core.hpp (C++) e de fcc220_reprodutibilidade.py.
!
! Compilar:  gfortran -O2 baumplane_mef.F90 -o baumplane_mef_f
! Precisão estendida: quádrupla (113 bits de mantissa) se o gfortran tiver; senão 80 bits (64 bits de mantissa):
!            gfortran -O2 -DESTENDIDA baumplane_mef.F90 -o baumplane_mef_f_q
! Rodar:     baumplane_mef_f                     (FCC220-2V, estagnação)
!            baumplane_mef_f saida=html          (também grava baumplane_relatorio.html com gráficos,
!                                                 arquivos .dat e o script baumplane.gp do gnuplot)
!            baumplane_mef_f caso=op             (operação com a placa a 60 °C)
!            baumplane_mef_f eg=4 gap=16 eiso=40 (troca parâmetros: chave=valor)
module baumplane
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
#ifdef ESTENDIDA
  ! quádrupla (33 dígitos) se o compilador tiver; senão a estendida de 80 bits (18 dígitos)
  integer, parameter :: wp = merge(selected_real_kind(33), selected_real_kind(18), selected_real_kind(33) > 0)
#else
  integer, parameter :: wp = real64
#endif
  real(wp), parameter :: SIG = 5.670374e-8_wp, PI = 3.14159265358979323846264338327950_wp
  integer, parameter :: NMAX = 64, NODES = 13*NMAX + 1, NPAR = 16
  ! parâmetros (espessuras em mm), na mesma ordem do C++ e do app
  integer, parameter :: G = 1, TA = 2, BETA = 3, HW = 4, HI = 5, EG = 6, KG = 7, EPSG = 8, TAU = 9, &
                        GAP = 10, EABS = 11, KABS = 12, ALPHA = 13, EPS = 14, EISO = 15, KISO = 16
  character(len=5), parameter :: NOMES(NPAR) = [character(len=5) :: 'G', 'Ta', 'beta', 'hw', 'hi', 'eg', 'kg', &
                                  'epsg', 'tau', 'gap', 'eabs', 'kabs', 'alpha', 'eps', 'eiso', 'kiso']
  real(wp), parameter :: PADRAO(NPAR) = [800._wp, 25._wp, 45._wp, 10._wp, 5._wp, 3.2_wp, 1.0_wp, 0.88_wp, &
                                         0.91_wp, 25._wp, 0.4_wp, 237._wp, 0.95_wp, 0.05_wp, 25._wp, 0.040_wp]  ! FCC220-2V

  type :: analitica                 ! solução por regiões; c = A1 B1 A2 B2 A3 B3
    real(wp) :: k1, k2, k3, L1, L2, L3, he, hi, Tinf, q, c(6)
  end type
  type :: caso
    type(analitica) :: s
    real(wp) :: hg, qExtr, Sabs, Tp
  end type

contains

  ! ---------- correlações (Duffie e Beckman, 2013, cap. 3) ----------
  real(wp) function hGap(Tp, Tg, L, b)          ! Hollands, eq. 3.11.4
    real(wp), intent(in) :: Tp, Tg, L, b
    real(wp) :: TmK, k, nu, a, dT, Ra, c, t1, Nus   ! Fortran não diferencia maiúsculas: Nu seria nu
    TmK = (Tp + Tg)/2 + 273.15_wp
    k = 0.0241_wp*potr(TmK/273.15_wp, 0.81_wp)
    nu = 1.33e-5_wp*potr(TmK/273.15_wp, 1.75_wp)
    a = nu/0.71_wp
    dT = max(abs(Tp - Tg), 1e-6_wp)
    Ra = 9.81_wp/TmK*dT*L*L*L/(nu*a)
    c = cosr(b*PI/180)
    t1 = max(0._wp, 1 - 1708/(Ra*c))
    Nus = 1 + 1.44_wp*(1 - 1708*potr(sinr(1.8_wp*b*PI/180), 1.6_wp)/(Ra*c))*t1 + max(0._wp, potr(Ra*c/5830, 1._wp/3) - 1)
    hGap = Nus*k/L
  end function

  real(wp) function hRad(T1, T2, e1, e2)       ! placas paralelas, eq. 3.10.2
    real(wp), intent(in) :: T1, T2, e1, e2
    real(wp) :: a, b
    a = T1 + 273.15_wp; b = T2 + 273.15_wp
    hRad = SIG*(a*a + b*b)*(a + b)/(1/e1 + 1/e2 - 1)
  end function

  ! ---------- potência, seno, cosseno e raiz ----------
  ! Normalmente são as funções intrínsecas. Com -DMATEMATICA_PROPRIA usam séries em aritmética pura, para o
  ! núcleo em WebAssembly com precisão quádrupla (o flang 18 ainda não tem biblioteca matemática de kind 16).
  real(wp) function potr(x, y)                  ! x**y, x > 0
    real(wp), intent(in) :: x, y
#ifdef MATEMATICA_PROPRIA
    potr = mexp(y*mlog(x))
#else
    potr = x**y
#endif
  end function
  real(wp) function cosr(x)
    real(wp), intent(in) :: x
#ifdef MATEMATICA_PROPRIA
    cosr = taylor(x, 0)
#else
    cosr = cos(x)
#endif
  end function
  real(wp) function sinr(x)
    real(wp), intent(in) :: x
#ifdef MATEMATICA_PROPRIA
    sinr = taylor(x, 1)
#else
    sinr = sin(x)
#endif
  end function
  real(wp) function raiz(x)
    real(wp), intent(in) :: x
#ifdef MATEMATICA_PROPRIA
    integer :: i
    raiz = 0
    if (x <= 0) return
    raiz = real(sqrt(real(x, real64)), wp)    ! chute em dupla e Newton (cada passo dobra os dígitos)
    do i = 1, 3
      raiz = (raiz + x/raiz)/2
    end do
#else
    raiz = sqrt(x)
#endif
  end function
#ifdef MATEMATICA_PROPRIA
  real(wp) function taylor(x, m)                ! m = 0: cos x, m = 1: sen x; |x| até uns 2*pi
    real(wp), intent(in) :: x
    integer, intent(in) :: m
    real(wp) :: termo
    integer :: n
    termo = 1
    if (m == 1) termo = x
    taylor = termo
    do n = m + 2, 90, 2
      termo = -termo*x*x/(n*(n - 1))
      taylor = taylor + termo
    end do
  end function
  real(wp) function mlog(x)                     ! ln x = k ln 2 + 2 atanh((y - 1)/(y + 1)), y em [1, 2)
    real(wp), intent(in) :: x
    real(wp), parameter :: LN2 = 0.693147180559945309417232121458176568_wp
    real(wp) :: y, z, z2, termo
    integer :: k, n
    y = x; k = 0
    do while (y >= 2)
      y = y/2; k = k + 1
    end do
    do while (y < 1)
      y = y*2; k = k - 1
    end do
    z = (y - 1)/(y + 1); z2 = z*z; termo = z; mlog = 0
    do n = 1, 81, 2
      mlog = mlog + termo/n; termo = termo*z2
    end do
    mlog = k*LN2 + 2*mlog
  end function
  real(wp) function mexp(x)                     ! e**x = 2**k e**r, e**r pela série de r/16 e quatro quadrados
    real(wp), intent(in) :: x
    real(wp), parameter :: LN2 = 0.693147180559945309417232121458176568_wp
    real(wp) :: r, termo
    integer :: k, n
    k = int(x/LN2); r = (x - k*LN2)/16
    termo = 1; mexp = 1
    do n = 1, 30
      termo = termo*r/n; mexp = mexp + termo
    end do
    do n = 1, 4
      mexp = mexp*mexp
    end do
    do n = 1, abs(k)
      if (k > 0) then
        mexp = mexp*2
      else
        mexp = mexp/2
      end if
    end do
  end function
#endif

  ! ---------- solução analítica: sistema 6×6 (Apêndice A) ----------
  real(wp) function T_an(s, x)
    type(analitica), intent(in) :: s
    real(wp), intent(in) :: x
    if (x <= s%L1) then
      T_an = s%c(1)*x + s%c(2)
    else if (x <= s%L2) then
      T_an = -s%q/(2*s%k2)*x*x + s%c(3)*x + s%c(4)
    else
      T_an = s%c(5)*x + s%c(6)
    end if
  end function

  function analitica6(k1, k2, k3, L1, L2, L3, he, hi, Tinf, q) result(s)
    real(wp), intent(in) :: k1, k2, k3, L1, L2, L3, he, hi, Tinf, q
    type(analitica) :: s
    real(wp) :: M(6, 6), b(6), f, tb
    integer :: c, r, p, k
    M = 0                                       ! (sem reshape: o núcleo não depende da biblioteca do Fortran)
    M(1, 1) = k1; M(1, 2) = -he
    M(2, 5) = -(k3 + hi*L3); M(2, 6) = -hi
    M(3, 1) = L1; M(3, 2) = 1; M(3, 3) = -L1; M(3, 4) = -1
    M(4, 1) = -k1; M(4, 3) = k2
    M(5, 3) = L2; M(5, 4) = 1; M(5, 5) = -L2; M(5, 6) = -1
    M(6, 3) = -k2; M(6, 5) = k3
    b(1) = -he*Tinf; b(2) = -hi*Tinf; b(3) = -q*L1*L1/(2*k2); b(4) = q*L1; b(5) = q*L2*L2/(2*k2); b(6) = -q*L2
    do c = 1, 6                                 ! Gauss com pivotamento parcial
      p = c
      do r = c + 1, 6
        if (abs(M(r, c)) > abs(M(p, c))) p = r
      end do
      do k = 1, 6
        tb = M(c, k); M(c, k) = M(p, k); M(p, k) = tb
      end do
      tb = b(c); b(c) = b(p); b(p) = tb
      do r = c + 1, 6
        f = M(r, c)/M(c, c)
        do k = c, 6
          M(r, k) = M(r, k) - f*M(c, k)
        end do
        b(r) = b(r) - f*b(c)
      end do
    end do
    do r = 6, 1, -1
      tb = b(r)
      do k = r + 1, 6
        tb = tb - M(r, k)*s%c(k)
      end do
      s%c(r) = tb/M(r, r)
    end do
    s%k1 = k1; s%k2 = k2; s%k3 = k3; s%L1 = L1; s%L2 = L2; s%L3 = L3
    s%he = he; s%hi = hi; s%Tinf = Tinf; s%q = q
  end function

  ! ---------- MEF linear (Apêndice C): tridiagonal, algoritmo de Thomas ----------
  integer function malha(L1, L2, L3, n, x)      ! 4n / n / 8n elementos; devolve o nº de nós
    real(wp), intent(in) :: L1, L2, L3
    integer, intent(in) :: n
    real(wp), intent(out) :: x(:)
    integer :: i, N_
    N_ = 0
    do i = 0, 4*n;  N_ = N_ + 1; x(N_) = L1*i/(4*n);                end do
    do i = 1, n;    N_ = N_ + 1; x(N_) = L1 + (L2 - L1)*i/n;        end do
    do i = 1, 8*n;  N_ = N_ + 1; x(N_) = L2 + (L3 - L2)*i/(8*n);    end do
    malha = N_
  end function

  subroutine mef(s, x, N, T)
    type(analitica), intent(in) :: s
    real(wp), intent(in) :: x(:)
    integer, intent(in) :: N
    real(wp), intent(out) :: T(:)
    real(wp) :: a(NODES), d(NODES), c(NODES), F(NODES), xm, Le, k, qv, st, fe, m
    integer :: e, i
    a(1:N) = 0; d(1:N) = 0; c(1:N) = 0; F(1:N) = 0
    do e = 1, N - 1
      xm = (x(e) + x(e + 1))/2; Le = x(e + 1) - x(e)
      if (xm < s%L1) then; k = s%k1; else if (xm < s%L2) then; k = s%k2; else; k = s%k3; end if
      qv = 0; if (xm > s%L1 .and. xm < s%L2) qv = s%q
      st = k/Le; fe = qv*Le/2
      d(e) = d(e) + st; d(e + 1) = d(e + 1) + st; c(e) = c(e) - st; a(e + 1) = a(e + 1) - st
      F(e) = F(e) + fe; F(e + 1) = F(e + 1) + fe
    end do
    d(1) = d(1) + s%he; F(1) = F(1) + s%he*s%Tinf          ! Robin na frente
    d(N) = d(N) + s%hi; F(N) = F(N) + s%hi*s%Tinf          ! Robin no fundo
    do i = 2, N
      m = a(i)/d(i - 1); d(i) = d(i) - m*c(i - 1); F(i) = F(i) - m*F(i - 1)
    end do
    T(N) = F(N)/d(N)
    do i = N - 1, 1, -1
      T(i) = (F(i) - c(i)*T(i + 1))/d(i)
    end do
  end subroutine

  ! ---------- montagem com ponto fixo dos coeficientes ----------
  subroutine montar(P, qExtr, C)                ! (sub-rotina: sem cópia de tipo derivado pela biblioteca do Fortran)
    real(wp), intent(in) :: P(NPAR), qExtr
    type(caso), intent(out) :: C
    real(wp) :: eg_, Lg, ea, L1, L2, L3, Sabs, qv, Tp, Tg, To, hg, k1, A, B, he, nTp, nTo, nTg, dd
    type(analitica) :: s
    integer :: it
    eg_ = P(EG)/1000; Lg = P(GAP)/1000; ea = P(EABS)/1000
    L1 = eg_ + Lg; L2 = L1 + ea; L3 = L2 + P(EISO)/1000
    Sabs = P(G)*1.01_wp*P(TAU)*P(ALPHA); qv = (Sabs - qExtr)/ea
    Tp = 80; Tg = 40; To = 35
    do it = 1, 200
      hg = hGap(Tp, Tg, Lg, P(BETA)) + hRad(Tp, Tg, P(EPS), P(EPSG))
      k1 = L1/(eg_/P(KG) + 1/hg)
      A = To + 273.15_wp; B = P(TA) + 273.15_wp
      he = P(HW) + P(EPSG)*SIG*(A*A + B*B)*(A + B)          ! Robin + radiação linearizada
      s = analitica6(k1, P(KABS), P(KISO), L1, L2, L3, he, P(HI), P(TA), qv)
      nTp = T_an(s, L1); nTo = T_an(s, 0._wp); nTg = nTo + abs(-s%k1*s%c(1))*eg_/P(KG)
      dd = max(abs(nTp - Tp), abs(nTg - Tg), abs(nTo - To))
      Tp = nTp; Tg = nTg; To = nTo
      if (dd < 1e-12_wp) exit
    end do
    C%s = s; C%hg = hg; C%qExtr = qExtr; C%Sabs = Sabs; C%Tp = Tp
  end subroutine
end module baumplane

#ifndef SO_MODULO
! (baumplane_wasm.F90 define SO_MODULO e aproveita só o módulo acima)
program baumplane_mef_f
  use baumplane
  implicit none
  real(wp) :: P(NPAR), lo, hi_, m, qf, qb, dT, nrm, nod, e2, Er, Eant, p_, h, gp(3), gw(3)
  real(wp) :: x(NODES), T(NODES)
  type(caso) :: C
  character(len=64) :: arg, k, v, prec
  integer :: i, j, eq, n, N_, e
  logical :: op, piso, html
  integer :: nr, rn(8), rNn(8)                          ! linhas da tabela de convergência, para o relatório
  real(wp) :: rh(8), rnod(8), rEr(8), rp(8)
  logical :: rpiso(8)
  integer :: u                                          ! relatório: unidade do arquivo e escalas dos gráficos
  real(wp) :: ymin, ymax, lh0, lh1
  real(wp), parameter :: FW = 640, FH = 300, FL = 56, FR = 16, FT = 14, FB = 44   ! figura (px)

  P = PADRAO; op = .false.; html = .false.; nr = 0
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    eq = index(arg, '=')
    if (eq == 0) then
      print '(a)', 'argumento invalido: '//trim(arg)//' (use chave=valor)'; stop 1
    end if
    k = arg(:eq - 1); v = arg(eq + 1:)
    if (k == 'caso') then
      op = (v == 'op'); cycle
    end if
    if (k == 'saida') then
      html = (v == 'html'); cycle
    end if
    j = 0                                        ! (findloc só existe a partir do gfortran 9)
    do eq = 1, NPAR
      if (NOMES(eq) == k) j = eq
    end do
    if (j == 0) then
      print '(a)', 'parametro desconhecido: '//trim(k); stop 1
    end if
    read (v, *) P(j)
  end do

  call montar(P, 0._wp, C)
  if (op) then                                   ! extração que deixa a placa a 60 °C (bissecção)
    if (C%Tp <= 60) then
      print '(a)', 'Em estagnacao a placa nao chega a 60 C; mostrando estagnacao.'
    else
      lo = 0; hi_ = C%Sabs
      do i = 1, 80
        m = (lo + hi_)/2; call montar(P, m, C)
        if (C%Tp > 60) then; lo = m; else; hi_ = m; end if
      end do
    end if
  end if

  associate (s => C%s)
    qf = s%k1*s%c(1); qb = -s%k3*s%c(5); dT = C%Tp - P(TA)
    write (prec, '(a,i0,a)') 'mantissa de ', digits(1._wp), ' bits'   ! 53 dupla, 64 estendida, 113 quadrupla
    if (op) then
      print '(a)', 'BaumPlane_colector: MEF 1D em Fortran (operacao, placa a 60 C, '//trim(prec)//')'
    else
      print '(a)', 'BaumPlane_colector: MEF 1D em Fortran (estagnacao, '//trim(prec)//')'
    end if
    print '(a,f0.3,a,f0.6,a)', '  placa T(L1)         = ', C%Tp, ' C   (Tmax no absorvedor ', &
          s%k2*s%c(3)**2/(2*s%q) + s%c(4), ' C)'
    print '(a,f0.3,a,f0.4,a,f0.2,a)', '  h camara            = ', C%hg, ' W/m2K, k1 equivalente = ', s%k1, &
          ' W/mK, h externo = ', s%he, ' W/m2K'
    print '(a,f0.2,a,f0.2,a,f0.2,a,f0.2,a,es8.1,a)', '  balanco             : frente ', qf, ' + fundo ', qb, &
          ' + fluido ', C%qExtr, ' = S ', C%Sabs, ' W/m2 (residuo ', abs(qf + qb + C%qExtr - C%Sabs), ')'
    print '(a,f0.4,a,f0.4,a,f0.4,a)', '  U frente + fundo    = ', qf/dT, ' + ', qb/dT, ' = ', (qf + qb)/dT, ' W/m2K'
    print '(a)', ''

    ! norma de T: Gauss 3 pontos é exato para T² (grau ≤ 4) em cada região
    gp = [-sqrt(0.6_wp), 0._wp, sqrt(0.6_wp)]; gw = [5._wp/9, 8._wp/9, 5._wp/9]
    nrm = sqrt(quadT2(0._wp, s%L1) + quadT2(s%L1, s%L2) + quadT2(s%L2, s%L3))

    print '(a)', '   n    nos   h_II (mm)   erro nodal max (K)   E_r (L2)      p'
    n = 1; Eant = 0; piso = .false.
    do while (n <= NMAX)
      N_ = malha(s%L1, s%L2, s%L3, n, x)
      call mef(s, x, N_, T)
      nod = 0; e2 = 0
      do i = 1, N_
        nod = max(nod, abs(T(i) - T_an(s, x(i))))
      end do
      do e = 1, N_ - 1
        e2 = e2 + quadErr(x(e), x(e + 1), T(e), T(e + 1))
      end do
      Er = sqrt(e2)/nrm; h = (s%L2 - s%L1)/n*1000
      if (n == 1) then
        print '(2x,i3,2x,i5,3x,f9.5,3x,es18.2,3x,es9.2,7x,a)', n, N_, h, nod, Er, '-'
      else
        p_ = log(Eant/Er)/log(2._wp)
        piso = piso .or. abs(p_ - 2) >= 0.15_wp
        if (piso) then
          print '(2x,i3,2x,i5,3x,f9.5,3x,es18.2,3x,es9.2,3x,f5.2,a)', n, N_, h, nod, Er, p_, '  (arredondamento)'
        else
          print '(2x,i3,2x,i5,3x,f9.5,3x,es18.2,3x,es9.2,3x,f5.2)', n, N_, h, nod, Er, p_
        end if
      end if
      nr = nr + 1; rn(nr) = n; rNn(nr) = N_; rh(nr) = h; rnod(nr) = nod; rEr(nr) = Er
      rp(nr) = 0; rpiso(nr) = piso; if (n > 1) rp(nr) = p_
      Eant = Er; n = 2*n
    end do
  end associate
  if (html) call relatorio()

contains
  ! ---------- relatório: HTML com gráficos SVG + arquivos .dat + script do gnuplot ----------
  subroutine relatorio()
    integer :: i, N1
    real(wp) :: x1(NODES), T1(NODES), xx, tmin, tmax
    character(len=*), parameter :: LARANJA = '#d95926', AZUL = '#2a78d6'
    associate (s => C%s)
      N1 = malha(s%L1, s%L2, s%L3, 1, x1)
      call mef(s, x1, N1, T1)
      ! --- dados para gnuplot ---
      open (newunit=u, file='baumplane_perfil.dat', status='replace')
      write (u, '(a)') '# x (mm)   T analitica (C)'
      do i = 0, 400
        xx = s%L3*i/400; write (u, '(f12.5,1x,f14.6)') xx*1000, T_an(s, xx)
      end do
      close (u)
      open (newunit=u, file='baumplane_nos.dat', status='replace')
      write (u, '(a)') '# x (mm)   T MEF n=1 (C)'
      do i = 1, N1
        write (u, '(f12.5,1x,f14.6)') x1(i)*1000, T1(i)
      end do
      close (u)
      open (newunit=u, file='baumplane_convergencia.dat', status='replace')
      write (u, '(a)') '# n   nos   h_II (mm)   erro nodal max (K)   E_r (L2)   p'
      do i = 1, nr
        write (u, '(i4,1x,i5,1x,f10.5,2(1x,es12.4),1x,f8.3)') rn(i), rNn(i), rh(i), rnod(i), rEr(i), rp(i)
      end do
      close (u)
      open (newunit=u, file='baumplane.gp', status='replace')
      write (u, '(a)') '# gnuplot baumplane.gp  ->  baumplane_perfil.png e baumplane_convergencia.png'
      write (u, '(a)') 'set terminal pngcairo size 900,500 font "Sans,12"'
      write (u, '(a)') 'set output "baumplane_perfil.png"'
      write (u, '(a)') 'set xlabel "x (mm)"; set ylabel "T (C)"; set grid; set key top left'
      write (u, '(a,f0.4,a,f0.4,a)') 'set arrow from ', s%L1*1000, ',graph 0 to ', s%L1*1000, ',graph 1 nohead dt 2'
      write (u, '(a,f0.4,a,f0.4,a)') 'set arrow from ', s%L2*1000, ',graph 0 to ', s%L2*1000, ',graph 1 nohead dt 2'
      write (u, '(a)') 'plot "baumplane_perfil.dat" w l lw 2 lc rgb "'//AZUL//'" t "Analitica (6x6)", \'
      write (u, '(a)') '     "baumplane_nos.dat" w p pt 7 ps 1.2 lc rgb "'//LARANJA//'" t "MEF, n = 1"'
      write (u, '(a)') 'unset arrow; set output "baumplane_convergencia.png"'
      write (u, '(a)') 'set logscale xy; set format y "10^{%L}"; set xlabel "h_{II} (mm)"; set ylabel "E_r (L2)"'
      write (u, '(a)') 'f(x) = A*x**2; A = 1; fit [0.1:0.4] f(x) "baumplane_convergencia.dat" u 3:5 via A'
      write (u, '(a)') 'plot "baumplane_convergencia.dat" u 3:5 w lp pt 7 lw 2 lc rgb "'//LARANJA//'" t "MEF", \'
      write (u, '(a)') '     f(x) w l dt 2 lc rgb "#666666" t "inclinacao 2 (h^2)"'
      close (u)

      ! --- HTML com SVG ---
      open (newunit=u, file='baumplane_relatorio.html', status='replace')
      write (u, '(a)') '<!doctype html><html lang="pt-BR"><head><meta charset="utf-8">'
      write (u, '(a)') '<meta name="viewport" content="width=device-width, initial-scale=1"><title>BaumPlane_MEF: relat'// &
            'orio Fortran</title>'
      write (u, '(a)') '<style>body{font:14px/1.5 system-ui,sans-serif;max-width:720px;margin:24px auto;padding:0 16px;'// &
            'color:#1d2421}'
      write (u, '(a)') 'h1{font:600 26px Georgia,serif;margin:0}table{border-collapse:collapse;width:100%;font:13px ui-'// &
            'monospace,monospace}'
      write (u, '(a)') 'td,th{padding:4px 8px;border-bottom:1px solid #ddd;text-align:right}th{color:#5d6863;font-weight:500}'
      write (u, '(a)') '.m{color:#5d6863}svg{max-width:100%;height:auto}svg text{font:11px system-ui;fill:#5d6863}</sty'// &
            'le></head><body>'
      write (u, '(a)') '<h1>BaumPlane_MEF</h1><p class="m">Relatorio gerado pelo programa Fortran baumplane_mef.F90 ('// &
                       trim(prec)//')</p>'
      write (u, '(a)') '<table><tr><th style="text-align:left">Grandeza</th><th>Valor</th></tr>'
      write (u, '(a,f0.3,a)') '<tr><td style="text-align:left">T placa</td><td>', C%Tp, ' &deg;C</td></tr>'
      write (u, '(a,f0.2,a,f0.2,a,f0.2,a)') '<tr><td style="text-align:left">frente + fundo + fluido</td><td>', &
            qf, ' + ', qb, ' + ', C%qExtr, ' W/m&sup2;</td></tr>'
      write (u, '(a,f0.2,a)') '<tr><td style="text-align:left">S absorvido</td><td>', C%Sabs, ' W/m&sup2;</td></tr>'
      write (u, '(a,f0.4,a)') '<tr><td style="text-align:left">U frente + fundo</td><td>', (qf + qb)/dT, ' W/(m&sup2;K)'// &
            '</td></tr>'
      write (u, '(a,f0.3,a,f0.4,a)') '<tr><td style="text-align:left">h camara / k1 equivalente</td><td>', C%hg, &
            ' / ', s%k1, '</td></tr></table>'

      ! perfil de temperatura
      tmin = minval([(T_an(s, s%L3*i/400), i = 0, 400)]); tmax = maxval([(T_an(s, s%L3*i/400), i = 0, 400)])
      ymin = 20*floor(tmin/20); ymax = 20*ceiling(tmax/20)
      write (u, '(a)') '<h2>Perfil de temperatura</h2><p class="m"><span style="color:'//AZUL// &
                       '">&#9473;</span> analitica (6x6) &nbsp; <span style="color:'//LARANJA//'">&#9679;</span> nos do'// &
                             ' MEF, n = 1</p>'
      write (u, '(a,i0,a,i0,a)') '<svg viewBox="0 0 ', nint(FW), ' ', nint(FH), '">'
      call regiao(0._wp, s%L1, 'I'); call regiao(s%L1, s%L2, 'II'); call regiao(s%L2, s%L3, 'III')
      xx = ymin
      do while (xx <= ymax + 1e-9_wp)
        write (u, '(a,f0.1,a,f0.1,a,f0.1,a,f0.1,a)') '<line x1="', FL, '" x2="', FW - FR, '" y1="', py(xx), '" y2="', py(xx), &
              '" stroke="#e3e7e4"/>'
        write (u, '(a,f0.1,a,f0.1,a,i0,a)') '<text x="', FL - 6, '" y="', py(xx) + 4, '" text-anchor="end">', nint(xx), &
              '</text>'
        xx = xx + 20
      end do
      write (u, '(a)', advance='no') '<path fill="none" stroke="'//AZUL//'" stroke-width="2" d="'
      do i = 0, 400
        xx = s%L3*i/400
        if (i == 0) then
          write (u, '(a,f0.1,a,f0.1)', advance='no') 'M', px(xx), ',', py(T_an(s, xx))
        else
          write (u, '(a,f0.1,a,f0.1)', advance='no') ' L', px(xx), ',', py(T_an(s, xx))
        end if
      end do
      write (u, '(a)') '"/>'
      do i = 1, N1
        write (u, '(a,f0.1,a,f0.1,a)') '<circle cx="', px(x1(i)), '" cy="', py(T1(i)), '" r="3.5" fill="'//LARANJA//'"/>'
      end do
      write (u, '(a,f0.1,a,f0.1,a,f0.1,a)') '<text x="', (FW + FL)/2, '" y="', FH - 8, '" text-anchor="middle">x (mm), es'// &
            'pessura total ', &
            s%L3*1000, ' mm</text></svg>'

      ! convergência em escala log-log
      lh0 = log10(rh(nr)) - 0.15_wp; lh1 = log10(rh(1)) + 0.15_wp
      ymax = ceiling(log10(maxval(rEr(1:nr))))         ! a faixa inclui o fim da reta de inclinação 2
      ymin = floor(min(log10(minval(rEr(1:nr))), log10(rEr(1)) - 2*(log10(rh(1)) - log10(rh(nr)))))
      write (u, '(a)') '<h2>Convergencia do MEF (log-log)</h2><p class="m"><span style="color:'//LARANJA// &
                       '">&#9679;</span> E<sub>r</sub> (L&sup2;) &nbsp; <span>- - -</span> inclinacao 2 (h&sup2;)</p>'
      write (u, '(a,i0,a,i0,a)') '<svg viewBox="0 0 ', nint(FW), ' ', nint(FH), '">'
      xx = ymin
      do while (xx <= ymax + 1e-9_wp)
        write (u, '(a,f0.1,a,f0.1,a,f0.1,a,f0.1,a)') '<line x1="', FL, '" x2="', FW - FR, '" y1="', ly(xx), '" y2="', ly(xx), &
              '" stroke="#e3e7e4"/>'
        write (u, '(a,f0.1,a,f0.1,a,i0,a)') '<text x="', FL - 6, '" y="', ly(xx) + 4, '" text-anchor="end">1e', &
              nint(xx), '</text>'
        xx = xx + 1
      end do
      write (u, '(a,f0.1,a,f0.1,a,f0.1,a,f0.1,a)') '<line x1="', lx(log10(rh(1))), '" y1="', ly(log10(rEr(1))), &
            '" x2="', lx(log10(rh(nr))), '" y2="', ly(log10(rEr(1)) - 2*(log10(rh(1)) - log10(rh(nr)))), &
            '" stroke="#666" stroke-dasharray="5 4"/>'
      write (u, '(a)', advance='no') '<path fill="none" stroke="'//LARANJA//'" stroke-width="2" d="'
      do i = 1, nr
        if (i == 1) then
          write (u, '(a,f0.1,a,f0.1)', advance='no') 'M', lx(log10(rh(i))), ',', ly(log10(rEr(i)))
        else
          write (u, '(a,f0.1,a,f0.1)', advance='no') ' L', lx(log10(rh(i))), ',', ly(log10(rEr(i)))
        end if
      end do
      write (u, '(a)') '"/>'
      do i = 1, nr
        write (u, '(a,f0.1,a,f0.1,a,f0.1,a,f0.1,a,f0.4,a)') '<circle cx="', lx(log10(rh(i))), '" cy="', ly(log10(rEr(i))), &
              '" r="4" fill="'//LARANJA//'"/><text x="', lx(log10(rh(i))), '" y="', FH - FB + 16, '" text-anchor="middle'// &
                    '">', rh(i), '</text>'
      end do
      write (u, '(a,f0.1,a,f0.1,a)') '<text x="', (FW + FL)/2, '" y="', FH - 8, '" text-anchor="middle">h_II (mm)</text></svg>'
      write (u, '(a)') '<table><tr><th>n</th><th>nos</th><th>h_II (mm)</th><th>erro nodal (K)</th><th>E_r (L2)</th><th>'// &
            'p</th></tr>'
      do i = 1, nr
        if (i == 1) then
          write (u, '(a,i0,a,i0,a,f0.5,a,es9.2,a,es9.2,a)') '<tr><td>', rn(i), '</td><td>', rNn(i), '</td><td>', rh(i), &
                '</td><td>', rnod(i), '</td><td>', rEr(i), '</td><td>-</td></tr>'
        else
          write (u, '(a,i0,a,i0,a,f0.5,a,es9.2,a,es9.2,a,f5.2,a)') '<tr><td>', rn(i), '</td><td>', rNn(i), '</td><td>', rh(i), &
                '</td><td>', rnod(i), '</td><td>', rEr(i), '</td><td>', rp(i), merge(' (arred.)', '         ', &
                      rpiso(i))//'</td></tr>'
        end if
      end do
      write (u, '(a)') '</table><p class="m">Arquivos para o gnuplot: baumplane_perfil.dat, baumplane_nos.dat, '// &
                       'baumplane_convergencia.dat e baumplane.gp.</p></body></html>'
      close (u)
    end associate
    print '(a)', ''
    print '(a)', '  gravados: baumplane_relatorio.html, baumplane_perfil.dat, baumplane_nos.dat,'
    print '(a)', '            baumplane_convergencia.dat, baumplane.gp (gnuplot baumplane.gp gera os PNG)'
  end subroutine

  ! escalas dos gráficos SVG do relatório
    real(wp) function px(xv)
      real(wp), intent(in) :: xv
      px = FL + xv/C%s%L3*(FW - FL - FR)
    end function
    real(wp) function py(tv)
      real(wp), intent(in) :: tv
      py = FT + (1 - (tv - ymin)/(ymax - ymin))*(FH - FT - FB)
    end function
    real(wp) function lx(lv)
      real(wp), intent(in) :: lv
      lx = FL + (lh1 - lv)/(lh1 - lh0)*(FW - FL - FR)          ! h decrescente para a direita
    end function
    real(wp) function ly(lv)
      real(wp), intent(in) :: lv
      ly = FT + (ymax - lv)/(ymax - ymin)*(FH - FT - FB)
    end function
    subroutine regiao(a, b, nome)
      real(wp), intent(in) :: a, b
      character(len=*), intent(in) :: nome
      write (u, '(a,f0.1,a,f0.1,a,f0.1,a,f0.1,a)') '<rect x="', px(a), '" y="', FT, '" width="', max(px(b) - px(a), 1._wp), &
            '" height="', FH - FT - FB, '" fill="#eef1ef"/>'
      write (u, '(a,f0.1,a,f0.1,a)') '<text x="', (px(a) + px(b))/2, '" y="', FT + 14, '" text-anchor="middle"'// &
            '>'//nome//'</text>'
    end subroutine

  real(wp) function quadT2(a, b)
    real(wp), intent(in) :: a, b
    integer :: jj
    real(wp) :: z, t
    quadT2 = 0
    do jj = 1, 3
      z = (b - a)/2*gp(jj) + (a + b)/2; t = T_an(C%s, z)
      quadT2 = quadT2 + gw(jj)*t*t
    end do
    quadT2 = (b - a)/2*quadT2
  end function
  real(wp) function quadErr(a, b, Ta_, Tb)
    real(wp), intent(in) :: a, b, Ta_, Tb
    integer :: jj
    real(wp) :: z, d
    quadErr = 0
    do jj = 1, 3
      z = (b - a)/2*gp(jj) + (a + b)/2
      d = Ta_ + (Tb - Ta_)*(z - a)/(b - a) - T_an(C%s, z)
      quadErr = quadErr + gw(jj)*d*d
    end do
    quadErr = (b - a)/2*quadErr
  end function
end program
#endif
