# BaumPlane_colector: MEF 1D em Fortran

`baumplane_mef.F90` reproduz em Fortran moderno (2008) o mesmo modelo e o mesmo algoritmo do núcleo C++ (`../cpp/baumplane_core.hpp`): o modelo por regiões, a solução analítica 6×6, o MEF linear com malha 4n/n/8n e o erro L² por quadratura de Gauss. É um arquivo só, sem dependências.

## Compilar

    gfortran -O2 baumplane_mef.F90 -o baumplane_mef_f
    gfortran -O2 -DESTENDIDA baumplane_mef.F90 -o baumplane_mef_f_q

O `-DESTENDIDA` usa a precisão quádrupla (113 bits de mantissa) quando o gfortran a oferece; se não oferecer, cai na estendida de 80 bits (64 bits de mantissa). A primeira linha da saída diz qual foi usada. A extensão `.F90`, com F maiúsculo, faz o gfortran passar o arquivo pelo pré-processador.

No Windows, o pacote WinLibs (`winget install BrechtSanders.WinLibs.POSIX.UCRT`) já traz o gfortran junto com o g++.

## Rodar

    baumplane_mef_f                       (FCC220-2V em estagnação)
    baumplane_mef_f caso=op               (operação com a placa a 60 °C)
    baumplane_mef_f eg=4 gap=16 eiso=40   (troca parâmetros: chave=valor)

## Gráficos

    baumplane_mef_f saida=html

Além da tabela no terminal, o programa grava na pasta atual:
- `baumplane_relatorio.html`: abre no navegador, com o perfil de temperatura (analítica e nós do MEF), a convergência em escala log-log com a reta de inclinação 2 e as tabelas. É SVG escrito pelo próprio Fortran, sem nenhuma biblioteca.
- `baumplane_perfil.dat`, `baumplane_nos.dat` e `baumplane_convergencia.dat`: colunas de texto para qualquer programa de gráficos.
- `baumplane.gp`: script do gnuplot. Rodar `gnuplot baumplane.gp` gera `baumplane_perfil.png` e `baumplane_convergencia.png`, e ajusta E = A·h² nas malhas grossas.

A pasta `exemplos/` traz essas saídas em dupla e em quádrupla precisão. Lado a lado, as duas mostram a curva de convergência se desviando da reta h² em dupla precisão e seguindo reta até n = 64 em quádrupla.

## Conferência
- **Dupla precisão:** `saida_double.txt` e `saida_operacao.txt` coincidem dígito a dígito com o C++, inclusive nos erros da faixa de arredondamento. As únicas diferenças são de formato do Fortran: `.1308` em vez de `0.1308` e `E-09` em vez de `e-09`.
- **Precisão quádrupla:** em `saida_quadrupla.txt`, o MEF tem p = 2,00 em todas as malhas até n = 64, e o erro nodal fica em 10⁻²⁷ a 10⁻²⁴ K. Nessa precisão o piso de arredondamento simplesmente não aparece, o que confirma que os desvios das malhas finas em dupla precisão são numéricos, não do método.

## WebAssembly (o app BaumPlane_MEF)

O app na web roda este mesmo módulo Fortran, compilado para WebAssembly (`baumplane_mef_f.wasm`, cerca de 9 kB). `baumplane_wasm.F90` inclui `baumplane_mef.F90` sem o programa principal (`SO_MODULO`) e acrescenta quatro funções `bind(c)` que o JavaScript chama: `bp_set(i, v)` põe o parâmetro i (1 a 16), `bp_run(modo, Tref, nmax)` resolve (modo 0 estagnação, 1 operação), `bp_get(i)` lê a saída i e `bp_T(x)` dá a solução analítica. O núcleo foi escrito sem nada que dependa da biblioteca de execução do Fortran (sem `reshape`, construtores de vetor ou cópia de tipo derivado), e por isso basta o `shim.c` com `memset` e `memcpy`. A receita de compilação está no cabeçalho de `baumplane_wasm.F90`. Os números do app conferem dígito a dígito com `saida_double.txt` e `saida_operacao.txt`.

O app tem também a precisão quádrupla (`baumplane_mef_fq.wasm`, cerca de 43 kB), compilada com `-DESTENDIDA -DMATEMATICA_PROPRIA`. O WebAssembly não tem `real(16)` em hardware, então a aritmética de 113 bits é feita em software pelas rotinas do compiler-rt, e potência, seno, cosseno e raiz vêm de séries escritas no próprio módulo (`potr`, `sinr`, `cosr`, `raiz`; sem a opção, são as intrínsecas). Com ela o app mostra p = 2,00 em todas as malhas até n = 64, como `saida_quadrupla.txt`; Tp, fluxos e E_r conferem com o `baumplane_mef_f_q`, e só as colunas no nível do arredondamento (erro nodal perto de 1e-27, resíduo do balanço) mudam no último dígito, porque as séries não são idênticas à libquadmath do gfortran.
