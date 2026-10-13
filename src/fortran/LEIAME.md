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

## Conferência
- **Dupla precisão:** `saida_double.txt` e `saida_operacao.txt` coincidem dígito a dígito com o C++, inclusive nos erros da faixa de arredondamento. As únicas diferenças são de formato do Fortran: `.1308` em vez de `0.1308` e `E-09` em vez de `e-09`.
- **Precisão quádrupla:** em `saida_quadrupla.txt`, o MEF tem p = 2,00 em todas as malhas até n = 64, e o erro nodal fica em 10⁻²⁷ a 10⁻²⁴ K. Nessa precisão o piso de arredondamento simplesmente não aparece, o que confirma que os desvios das malhas finas em dupla precisão são numéricos, não do método.
