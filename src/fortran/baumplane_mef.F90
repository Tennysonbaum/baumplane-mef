! BaumPlane_colector: MEF 1D em Fortran (modelo por regiões, Apêndices A e C do TCC).
! Mesma física e mesmo algoritmo de baumplane_core.hpp (C++) e de fcc220_reprodutibilidade.py.
!
! Compilar:  gfortran -O2 baumplane_mef.F90 -o baumplane_mef_f
! Precisão estendida: quádrupla (113 bits de mantissa) se o gfortran tiver; senão 80 bits (64 bits de mantissa):
!            gfortran -O2 -DESTENDIDA baumplane_mef.F90 -o baumplane_mef_f_q
! Rodar:     baumplane_mef_f                     (FCC220-2V, estagnação)
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
    k = 0.0241_wp*(TmK/273.15_wp)**0.81_wp
    nu = 1.33e-5_wp*(TmK/273.15_wp)**1.75_wp
    a = nu/0.71_wp
    dT = max(abs(Tp - Tg), 1e-6_wp)
    Ra = 9.81_wp/TmK*dT*L*L*L/(nu*a)
    c = cos(b*PI/180)
    t1 = max(0._wp, 1 - 1708/(Ra*c))
    Nus = 1 + 1.44_wp*(1 - 1708*sin(1.8_wp*b*PI/180)**1.6_wp/(Ra*c))*t1 + max(0._wp, (Ra*c/5830)**(1._wp/3) - 1)
    hGap = Nus*k/L
  end function

  real(wp) function hRad(T1, T2, e1, e2)       ! placas paralelas, eq. 3.10.2
    real(wp), intent(in) :: T1, T2, e1, e2
    real(wp) :: a, b
    a = T1 + 273.15_wp; b = T2 + 273.15_wp
    hRad = SIG*(a*a + b*b)*(a + b)/(1/e1 + 1/e2 - 1)
  end function

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
    real(wp) :: M(6, 6), b(6), f, t(6), tb
    integer :: c, r, p, k
    M = reshape([real(wp) :: k1, -he, 0, 0, 0, 0, &
                 0, 0, 0, 0, -(k3 + hi*L3), -hi, &
                 L1, 1, -L1, -1, 0, 0, &
                 -k1, 0, k2, 0, 0, 0, &
                 0, 0, L2, 1, -L2, -1, &
                 0, 0, -k2, 0, k3, 0], [6, 6], order=[2, 1])
    b = [-he*Tinf, -hi*Tinf, -q*L1*L1/(2*k2), q*L1, q*L2*L2/(2*k2), -q*L2]
    do c = 1, 6                                 ! Gauss com pivotamento parcial
      p = c
      do r = c + 1, 6
        if (abs(M(r, c)) > abs(M(p, c))) p = r
      end do
      t = M(c, :); M(c, :) = M(p, :); M(p, :) = t
      tb = b(c); b(c) = b(p); b(p) = tb
      do r = c + 1, 6
        f = M(r, c)/M(c, c)
        M(r, c:) = M(r, c:) - f*M(c, c:)
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
    real(wp) :: a(N), d(N), c(N), F(N), xm, Le, k, qv, st, fe, m
    integer :: e, i
    a = 0; d = 0; c = 0; F = 0
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
  function montar(P, qExtr) result(C)
    real(wp), intent(in) :: P(NPAR), qExtr
    type(caso) :: C
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
  end function
end module baumplane

program baumplane_mef_f
  use baumplane
  implicit none
  real(wp) :: P(NPAR), lo, hi_, m, qf, qb, dT, nrm, nod, e2, Er, Eant, p_, h, gp(3), gw(3)
  real(wp) :: x(NODES), T(NODES)
  type(caso) :: C
  character(len=64) :: arg, k, v, prec
  integer :: i, j, eq, n, N_, e
  logical :: op, piso

  P = PADRAO; op = .false.
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
    j = 0                                        ! (findloc só existe a partir do gfortran 9)
    do eq = 1, NPAR
      if (NOMES(eq) == k) j = eq
    end do
    if (j == 0) then
      print '(a)', 'parametro desconhecido: '//trim(k); stop 1
    end if
    read (v, *) P(j)
  end do

  C = montar(P, 0._wp)
  if (op) then                                   ! extração que deixa a placa a 60 °C (bissecção)
    if (C%Tp <= 60) then
      print '(a)', 'Em estagnacao a placa nao chega a 60 C; mostrando estagnacao.'
    else
      lo = 0; hi_ = C%Sabs
      do i = 1, 80
        m = (lo + hi_)/2; C = montar(P, m)
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
      Eant = Er; n = 2*n
    end do
  end associate

contains
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
