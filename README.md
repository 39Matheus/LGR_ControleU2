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
