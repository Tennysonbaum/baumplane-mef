# BaumPlane_colector: núcleo MEF em C++

O núcleo de cálculo fica em `baumplane_core.hpp`. Ele tem a mesma física do app e do `fcc220_reprodutibilidade.py`: o modelo por regiões, a solução analítica 6×6, o MEF linear com malha 4n/n/8n e o erro L² por quadratura de Gauss. O núcleo não usa nem a biblioteca padrão, e por isso o mesmo código compila de duas formas:

- `baumplane_mef.cpp`: programa de linha de comando (g++).
- `baumplane_wasm.cpp`: módulo WebAssembly, que é o que roda no app **BaumPlane_MEF**.

Guarde os três arquivos na mesma pasta.

## Compilar
No Windows, instale o g++ uma vez pelo cmd:

    winget install BrechtSanders.WinLibs.POSIX.UCRT

Feche e abra o cmd de novo. Depois:

    g++ -O2 -std=c++17 baumplane_mef.cpp -o baumplane_mef.exe
    g++ -O2 -std=c++17 -DESTENDIDA baumplane_mef.cpp -o baumplane_mef_ld.exe

No Linux ou Mac, os comandos são os mesmos, sem o `.exe`.

## Rodar

    baumplane_mef                       (FCC220-2V em estagnação)
    baumplane_mef caso=op               (operação com a placa a 60 °C)
    baumplane_mef eg=4 gap=16 eiso=40   (troca parâmetros: chave=valor)
    baumplane_mef_ld                    (precisão estendida, long double)

Parâmetros, com as unidades do app (espessuras em mm): G, Ta, beta, hw, hi, eg, kg, epsg, tau, gap, eabs, kabs, alpha, eps, eiso, kiso.

## Conferência
- `saida_double.txt` e `saida_operacao.txt` reproduzem o relatório do FCC220-2V: 522,8 + 175,7 W/m², U = 4,82 W/(m²·K) e p = 2,00.
- `saida_long_double.txt` mostra a precisão estendida empurrando o piso de arredondamento de n = 8 para n = 64, que é o mesmo recurso usado no TCC.

## WebAssembly (o app BaumPlane_MEF)
Basta o clang com o wasm-ld (pacotes `clang` e `lld` no Linux, ou o LLVM no Windows). Não precisa do Emscripten:

    clang++ --target=wasm32 -O2 -std=c++17 -nostdlib -fno-exceptions -fno-rtti -fno-threadsafe-statics -fno-builtin -Wl,--no-entry -Wl,--allow-undefined -Wl,--export-dynamic baumplane_wasm.cpp -o baumplane_mef.wasm

O módulo tem cerca de 11 kB. As funções pow, cos, sin e cbrt vêm do Math do JavaScript, e sqrt e fabs são instruções do próprio wasm. Em estagnação, a tabela de convergência no navegador sai dígito a dígito igual à do g++ nativo.

## Próximo passo
O MEF 2D da aleta e da borda entra no mesmo núcleo, e o app passa a usá-lo sem duplicar código.
