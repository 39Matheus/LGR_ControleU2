# LGR_ControleU2

Calculadora Streamlit para a 2ª unidade de **Projeto de Sistemas de Controle**.

## Escopo

A versão acompanha a apostila até o capítulo 4, **Projeto de Controladores pelo
Método do LGR**. O capítulo 5, *Aproximação Discreta de Funções de Transferência
Contínuas*, permanece fora do escopo.

A entrada é única e possui uma seleção de operação:

- **Projeto de controlador**: PD, PI e PID com zeros reais iguais.
- **Análise do LGR**: roteiro clássico em 12 passos.

A planta e a realimentação permanecem preenchidas ao alternar entre as
operações. Os quatro exercícios da primeira lista da Unidade 2 são **presets de
preenchimento**: eles apenas carregam os campos normais, que continuam
editáveis antes da execução.

## Método de resolução

A saída de projeto foi organizada para reprodução manual:

1. converter Mp/ts, ξ/ωn ou polos dados em polos dominantes;
2. calcular as contribuições angulares por relações trigonométricas;
3. obter a posição do(s) zero(s) pela condição de ângulo, exibindo
   `Σφ_zeros - Σθ_polos = (2q+1)180°` termo a termo;
4. aplicar a condição de módulo no formato da resolução manual: distâncias
   `Ai` aos polos, `Bi` aos zeros, `KT=ΠAi/ΠBi` e
   `KC=KT/|KG KH|`;
5. converter para Kp/Ki/Kd;
6. conferir polos de malha fechada e resposta temporal;
7. opcionalmente, em especificações por desigualdade, aplicar ajuste fino após
   o projeto de fronteira. Esse recurso fica desativado por padrão para manter
   a resolução alinhada ao cálculo manual de prova.

Internamente atan2 é usado somente para escolher corretamente o quadrante; a
apresentação ao usuário segue o desenvolvimento trigonométrico da apostila.

Antes da análise, G(s)H(s) é reduzida simbolicamente. Fatores comuns exatamente
cancelados são removidos, e polos/zeros repetidos são calculados primeiro por
álgebra exata, com fallback numérico. Isso evita falhas de convergência do
SymPy em casos como a Questão 3.

## Organização

~~~text
core/
  lgr.py          # análise matemática do LGR
  specs.py        # Mp, ts, ξ, ωn e polos desejados
  geometry.py     # geometria/trigonometria
  controllers.py  # projeto PD/PI/PID
  simulation.py   # malha fechada e verificação
  presets.py      # exercícios da lista
ui/
  lgr_view.py
  controller_view.py
  renderers.py
reports/
  pdf.py
tests/
~~~

calculations.py foi mantido apenas como camada de compatibilidade.

## Executar

~~~bash
pip install -r requirements.txt
streamlit run app.py
~~~

## Testar

~~~bash
pip install pytest
pytest -q
~~~


## Exportação

A resolução de projeto possui duas saídas:

- **Modo Prova (PDF)**: versão compacta para transcrição manual.
- **Resolução completa**: pode ser baixada como **HTML interativo**, com os
  gráficos Plotly navegáveis, ou como **PDF completo**, mantendo a mesma
  sequência de resolução com gráficos estáticos. O nome do arquivo pode ser
  definido antes do download.
