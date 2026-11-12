# LGR_ControleU2

Calculadora Streamlit para a 2ª unidade de **Projeto de Sistemas de Controle**.

## Escopo

A versão acompanha agora o conteúdo utilizado nas duas listas da Unidade 2:

- análise do **Lugar Geométrico das Raízes**;
- projeto de controladores pelo LGR;
- **aproximação discreta de funções de transferência contínuas** (capítulo 5).

A entrada é única e possui três operações:

- **Projeto de controlador**;
- **Análise do LGR (12 passos)**;
- **Discretização**.

A planta e a realimentação permanecem preenchidas ao alternar entre as operações.

## Controladores

O projeto não fica mais restrito a PD, PI e PID com zeros iguais.

Estruturas prontas:

- P;
- PD;
- PI;
- PID com zeros reais **iguais ou distintos**;
- controlador **Genérico** da forma `Gc(s) = Kc Π(s + zi) / Π(s + pj)`.

No modo genérico podem ser informados zeros e polos reais fixos e escolhida uma singularidade real incógnita (zero ou polo), calculada pela condição de ângulo. O ganho é calculado pela condição de módulo.

Isso permite, por exemplo, resolver diretamente controladores do tipo `Gc(s)=a/(s+b)`.

Para PID com zeros distintos, um zero é informado pelo usuário e o segundo é determinado pela condição de ângulo. Isso evita assumir `z1=z2` quando o enunciado não impõe essa restrição.

## Método de projeto

A saída foi organizada para reprodução manual:

1. converter Mp/ts, ξ/ωn ou polos dados em polos dominantes;
2. calcular as contribuições angulares por relações trigonométricas;
3. montar `Σφ_zeros - Σθ_polos = (2q+1)180°`;
4. calcular a singularidade incógnita do controlador, quando existir;
5. aplicar a condição de módulo com distâncias `Ai` e `Bi`;
6. obter `KT=ΠAi/ΠBi` e `KC=KT/|KG KH|`;
7. converter para Kp/Ki/Kd quando a estrutura for P/PI/PD/PID;
8. conferir polos de malha fechada e resposta temporal;
9. opcionalmente aplicar ajuste fino em especificações por desigualdade.

O ajuste fino continua **desativado por padrão** para preservar os valores de fronteira usados na resolução manual de prova.

## Discretização

São implementados os métodos do capítulo 5:

- **Euler / Forward**: `s ≈ (z-1)/T`
- **Backward**: `s ≈ (z-1)/(Tz)`
- **Tustin / Bilinear**: `s ≈ (2/T)(z-1)/(z+1)`
- **Invariância ao degrau**, usando equivalência com segurador de ordem zero.

A resolução mostra a função contínua, a relação de aproximação, a função discretizada, coeficientes, polos/zeros em `z` e a equação de diferenças.

A discretização pode ser usada de forma independente ou encadeada ao projeto de controlador. No projeto é possível escolher como alvo `Gc(s)`, `G(s)`, `G(s)H(s)` ou `G(s)H(s)Gc(s)`.

## Presets

Os presets apenas **preenchem os campos normais**, que continuam editáveis.

- quatro questões da **Lista U2.1**;
- **Lista U2.2 — Questão 1**: PID com `z1=z2`, seguido de Tustin com `T=2 s`;
- **Lista U2.2 — Questão 2**: controlador genérico `a/(s+b)`, seguido de discretização de `G(s)H(s)Gc(s)` por Euler com `T=1 s`.

## Organização

```text
core/
  lgr.py
  specs.py
  geometry.py
  controllers.py
  discretization.py
  simulation.py
  presets.py
ui/
  main_view.py
  lgr_view.py
  controller_view.py
  discretization_view.py
  renderers.py
reports/
  pdf.py
  html.py
  full_pdf.py
tests/
```

`calculations.py` permanece como camada de compatibilidade.

## Exportação

- **Modo Prova (PDF)**: versão compacta;
- **HTML completo**: gráficos Plotly interativos;
- **PDF completo**: mesma sequência de resolução com gráficos estáticos.

Quando a discretização está habilitada após o projeto, ela também aparece nos relatórios.

## Enquadramento dos gráficos

O padrão é **automático proporcional**: +20% sobre o maior |Re| em `x` e +10% sobre o maior |Im| em `y`, mantendo piso de ±10. Também existem ajustes proporcional e linear personalizados.

## Execução

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Testes

```bash
pip install pytest
pytest -q
```


## Unidade 3 — Lista 1: espaço de estados

O menu inicial agora permite alternar entre o módulo preservado da **Unidade 2** e o módulo da **Unidade 3 / Lista 1**.

A U3L1 possui presets editáveis para **todos os itens**:

- **Q1, sistemas 1 e 2**: polinômio característico, autovalores/estabilidade, matrizes de controlabilidade `U=[B AB ... A^(n−1)B]` e observabilidade `V=[C;CA;...;CA^(n−1)]`, postos e conclusões;
- **Q2 e Q3**: projeto de realimentação de estados `u(t)=Kx(t)` e projeto do observador usando as fórmulas de Ackermann, com `q(A)`, matriz inversa, ganhos exatos e conferência por polinômios e polos;
- **Q4**: realização de `G(s)=(3s²+4s−2)/(s³+3s²+7s+5)` na **forma canônica controlável** da apostila, projeto do observador `−1,−1,−2` e realimentação `−4,−4,−5`.

**Convenção do professor:** `u=Kx`, portanto `A_mf=A+BK` e `K=−[0 ... 1]U⁻¹q(A)`. O observador usa `A−LC`, com `L=q(A)V⁻¹[0 ... 1]ᵀ`. As equações seguem as seções 7.3–7.8 da apostila.

A interface aceita matrizes e coeficientes editáveis, polos reais e pares conjugados, e disponibiliza exportação da resolução em Markdown. Os cálculos utilizam **SymPy exato**, preservando frações. Sistemas não controláveis/observáveis recebem aviso explicativo.

Testes da nova lista: `pytest -q tests/test_state_space.py`.


### Resolução reproduzível da U3L1

A tela foi refinada para **fórmula → cálculo intermediário → resultado → conferência**, mantendo as matrizes **A, B e C visíveis apenas uma vez**, na conferência das entradas. Na Q1, mostra `sI−A`, desenvolvimento do determinante, produtos `AB`, `A²B`, `CA`, `CA²`, `U`, `V`, determinantes e postos. Nas questões de Ackermann, mostra o polinômio desejado fatorado e expandido, `q(A)` numérico, a matriz de controlabilidade/observabilidade, sua inversa, o produto intermediário `eₙᵀU⁻¹` ou `V⁻¹eₙ`, o ganho `K` ou `L` e a conferência de `A+BK` ou `A−LC`.

Na forma canônica controlável, explicita os coeficientes `α` e `β` antes das matrizes. O botão de exportação Markdown reproduz esses mesmos passos; os dados alterados invalidam a solução anterior. Contas complementares são recolhidas para evitar poluição visual.


## Unidade 3, Lista 1 — relatórios completos por questão

A nova modalidade padrão **Questão completa (recomendado)** tem quatro presets, um para cada questão da U3L1:

1. **Q1:** a), b), c) para os **dois sistemas**, organizados em duas seções num único documento.
2. **Q2:** sistema 1 com a) realimentação de estados e b) observador de estados.
3. **Q3:** sistema 2 com c) realimentação de estados e d) observador de estados.
4. **Q4:** função de transferência da lista com a) forma canônica controlável, b) observador e c) realimentação.

Os documentos reproduzem os cálculos essenciais (fórmula, operações intermediárias, resultado, conferência), seguindo o capítulo 7 da apostila: **u=Kx => A+BK**, observador **A−LC**. As matrizes da entrada são exibidas uma só vez por sistema. Cada questão pode ser exportada para **PDF** ou **HTML** com MathJax; é possível também preparar um **ZIP com 4 PDFs e 4 HTMLs**. O modo **Item individual (editável)** preserva campos de matrizes, função de transferência, polos e seleção da operação, com exportação PDF/HTML específica. O PDF usa ReportLab e renderização matemática local, sem depender de serviços remotos nem LaTeX instalado no servidor; o HTML usa MathJax via CDN para fórmulas.

Execute `pytest -q tests/test_state_space_questions.py` para os testes de agrupamento, relatórios e ZIP.

## Unidade 3 — Lista 2: seguidores de referência (degrau)

A **U3L1 permanece integralmente disponível** no menu inicial. A nova opção **Unidade 3 — Lista 2** resolve as quatro questões da lista com a metodologia da apostila, seção 7.10.

- **Q1**: realização canônica controlável de `G(s)=(3s²+4s−2)/(s³+3s²+7s+5)` e seguidor de degrau com polos `−1±0,1i`, `−5`, `−12`.
- **Q2**: projeto no sistema `A=[[-6,1,0],[-11,0,1],[-6,0,0]]`, `B=[2,6,2]ᵀ`, `C=[1,0,0]`, com os mesmos quatro polos.
- **Q3**: realização canônica controlável de `G(s)=(2s+1)/(s(s+3))` e seguidor com os **três** polos `−1,−1,−1`.
- **Q4**: realização canônica **observável** de `G(s)=3/(s²+2s−4)` e seguidor com polos `−2,−2,−3`.

**Convenções idênticas à apostila** (eqs. 7.23–7.27): `e=y−r`, `z=dx/dt`, `w=du/dt`; matrizes aumentadas `A_a=[[0,C],[0,A]]`, `B_a=[[0],[B]]`; vetor de ganhos `K_a=−[0 ... 1]U_a⁻¹ q(A_a)=[k1,k2]`; controle `u(t)=k1∫(y−r)dt+k2x`; sistema fechado `[[A+Bk2,Bk1],[C,0]]`. O módulo verifica posto da controlabilidade aumentada, polos exatos e ganho estacionário unitário em degrau.

Cada um dos **4 presets de questão completa** apresenta dados uma única vez, forma canônica (quando solicitada), construção do sistema aumentado, produtos que formam a controlabilidade, polinômio desejado, `q(A_a)`, inversa, ganho, lei de controle e verificação. Os passos são idênticos na tela, PDF e HTML. O modo **Item individual (editável)** permite alterar os coeficientes da FT, matrizes e polos; gera PDF e HTML específicos. O pacote `.zip` contém exatamente **4 PDFs + 4 HTMLs**, nomeados `U3L2_Questao_1` a `U3L2_Questao_4`.

O motor PDF/HTML da U3L1 foi generalizado apenas para ler `list_label`, com padrão `U3L1`, de forma que os relatórios antigos e o ZIP antigo continuam com os nomes e conteúdos da U3L1.

Teste: `pytest -q tests/test_state_space_u3l2.py` ou `pytest -q` para a suíte completa.