"""Núcleo matemático do Lugar Geométrico das Raízes."""
import math
import warnings
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

try:
    import control as ctrl
except ImportError:
    ctrl = None

import numpy as np
import sympy as sp
import plotly.graph_objects as go

warnings.filterwarnings("ignore")


@dataclass
class VetorAnalise:
    origem: complex
    destino: complex
    dx: float
    dy: float
    magnitude: float
    angulo_deg: float


class AnalisadorLGR:
    """Núcleo de cálculos do Lugar Geométrico das Raízes.

    Esta classe não depende do Streamlit. Ela concentra os cálculos simbólicos,
    numéricos, tabelas e dados necessários para a interface apresentar uma
    resolução passo a passo.
    """

    def __init__(self, numG, denG, numH, denH):
        self.numG = np.asarray(numG, dtype=float)
        self.denG = np.asarray(denG, dtype=float)
        self.numH = np.asarray(numH, dtype=float)
        self.denH = np.asarray(denH, dtype=float)

        if ctrl is not None:
            self.G = ctrl.TransferFunction(numG, denG)
            self.H = ctrl.TransferFunction(numH, denH)
            self.GH = self.G * self.H
        else:
            self.G = None
            self.H = None
            self.GH = None

        self.s, self.K = sp.symbols("s K", real=True)
        self.w = sp.symbols("w", real=True)

        # Retira os coeficientes da função de transferência G(s)H(s) e os leva
        # para expressões simbólicas. K fica explicitamente separado.
        if ctrl is not None:
            num_arr = np.asarray(self.GH.num[0][0], dtype=float)
            den_arr = np.asarray(self.GH.den[0][0], dtype=float)
        else:
            # Convolução de numeradores e denominadores para formar G(s)H(s).
            num_arr = np.convolve(self.numG, self.numH)
            den_arr = np.convolve(self.denG, self.denH)

        self.num_expr = sp.Poly.from_list(num_arr.tolist(), gens=self.s).as_expr()
        self.den_expr = sp.Poly.from_list(den_arr.tolist(), gens=self.s).as_expr()

        self.num_expr = sp.factor(self.num_expr)
        self.den_expr = sp.factor(self.den_expr)
        self.P_expr = sp.cancel(self.num_expr / self.den_expr)

        # Ganho constante de P(s) quando escrito na forma fatorada:
        # P(s) = C * prod(s-z_i) / prod(s-p_i).
        num_poly = sp.Poly(self.num_expr, self.s)
        den_poly = sp.Poly(self.den_expr, self.s)
        self.ganho_constante = sp.simplify(num_poly.LC() / den_poly.LC())
        ganho_num = float(sp.N(self.ganho_constante))
        if abs(ganho_num) < 1e-15:
            raise ValueError("O ganho constante de P(s) não pode ser zero.")
        self.fase_constante_deg = 0.0 if ganho_num > 0 else 180.0

        # Polinômio característico: D(s) + K N(s) = 0.
        self.char_expr = sp.expand(self.den_expr + self.K * self.num_expr)
        self.char_poly = sp.Poly(self.char_expr, self.s)
        self.char_derivative = sp.expand(sp.diff(self.char_expr, self.s))

        self.polos: List[complex] = []
        self.zeros: List[complex] = []
        self.np = 0
        self.nz = 0

        self.real_segments: List[Tuple[float, float]] = []
        self.assymptote_center: Optional[float] = None
        self.assymptote_angles: List[float] = []
        self.breakaway_candidates: List[Dict[str, Any]] = []
        self.crossings: List[Dict[str, float]] = []
        self.departure_arrival: List[Dict[str, Any]] = []
        self.routh_used_epsilon = False

        self.min_x = -5.0
        self.max_x = 5.0
        self.span = 10.0
        self.fig = go.Figure()
        self._last_rlist = None
        self._last_klist = None
        self._configure_base_figure()

    # ------------------------------------------------------------------
    # Utilidades matemáticas / apresentação
    # ------------------------------------------------------------------
    @staticmethod
    def _complex_close_real(z: complex, tol=1e-7) -> bool:
        return abs(z.imag) < tol

    @staticmethod
    def _fmt_complex(z: complex, digits: int = 4) -> str:
        z = complex(z)
        if abs(z.imag) < 10 ** (-digits):
            return f"{z.real:.{digits}f}"
        sign = "+" if z.imag >= 0 else "-"
        return f"{z.real:.{digits}f} {sign} {abs(z.imag):.{digits}f}j"

    @staticmethod
    def _fmt_number(x: Any, digits: int = 4) -> str:
        try:
            xf = float(x)
            if abs(xf - round(xf)) < 10 ** (-digits):
                return str(int(round(xf)))
            return f"{xf:.{digits}f}"
        except Exception:
            return str(x)

    @staticmethod
    def _angle_deg(vector: complex) -> float:
        return math.degrees(math.atan2(vector.imag, vector.real))

    @staticmethod
    def _normalize_angle_360(angle: float) -> float:
        value = angle % 360.0
        if abs(value - 360.0) < 1e-10:
            value = 0.0
        return value

    def _configure_base_figure(self):
        self.fig.update_layout(
            xaxis_title="Eixo Real (Re)",
            yaxis_title="Eixo Imaginário (Im)",
            showlegend=True,
            hovermode="closest",
            height=550,
            plot_bgcolor="rgba(240, 240, 240, 0.5)",
            margin=dict(l=20, r=20, t=30, b=20),
            yaxis=dict(scaleanchor="x", scaleratio=1),
        )
        self.fig.add_hline(y=0, line_width=1.5, line_color="black", opacity=0.5)
        self.fig.add_vline(x=0, line_width=1.5, line_color="black", opacity=0.5)

    def _group_points(self, points: List[complex], tol: float = 1e-7):
        """Agrupa pontos coincidentes (dentro de uma tolerância) e informa multiplicidade."""
        groups = []
        for point in points:
            for group in groups:
                if abs(point - group["point"]) < tol:
                    group["multiplicity"] += 1
                    # Mantém uma média para reduzir ruído numérico.
                    m = group["multiplicity"]
                    group["point"] = ((m - 1) * group["point"] + point) / m
                    break
            else:
                groups.append({"point": point, "multiplicity": 1})
        return groups

    def _add_polos_zeros(self):
        pole_groups = self._group_points(self.polos)
        zero_groups = self._group_points(self.zeros)

        if pole_groups:
            norm = [g for g in pole_groups if g["multiplicity"] == 1]
            repeated = [g for g in pole_groups if g["multiplicity"] > 1]

            if norm:
                self.fig.add_trace(
                    go.Scatter(
                        x=[g["point"].real for g in norm],
                        y=[g["point"].imag for g in norm],
                        mode="markers",
                        marker=dict(symbol="x", size=14, color="red", opacity=1.0, line=dict(width=2)),
                        customdata=[g["multiplicity"] for g in norm],
                        hovertemplate="Pólo: %{x:.5g} + j%{y:.5g}<br>Multiplicidade: %{customdata}<extra></extra>",
                        name="Pólos",
                    )
                )

            if repeated:
                self.fig.add_trace(
                    go.Scatter(
                        x=[g["point"].real for g in repeated],
                        y=[g["point"].imag for g in repeated],
                        mode="markers",
                        marker=dict(symbol="x", size=14, color="red", opacity=0.42, line=dict(width=2)),
                        customdata=[g["multiplicity"] for g in repeated],
                        hovertemplate="Pólo múltiplo: %{x:.5g} + j%{y:.5g}<br>Multiplicidade: %{customdata}<extra></extra>",
                        name="Pólos múltiplos",
                    )
                )

        if zero_groups:
            self.fig.add_trace(
                go.Scatter(
                    x=[g["point"].real for g in zero_groups],
                    y=[g["point"].imag for g in zero_groups],
                    mode="markers",
                    marker=dict(symbol="circle-open", size=13, color="blue", line=dict(width=2)),
                    customdata=[g["multiplicity"] for g in zero_groups],
                    hovertemplate="Zero: %{x:.5g} + j%{y:.5g}<br>Multiplicidade: %{customdata}<extra></extra>",
                    name="Zeros",
                )
            )

    # ------------------------------------------------------------------
    # Passos 1 e 2
    # ------------------------------------------------------------------
    def calcular_polos_zeros(self) -> Dict[str, Any]:
        raizes_num = sp.roots(self.num_expr, self.s)
        raizes_den = sp.roots(self.den_expr, self.s)

        self.zeros = [complex(z.evalf()) for z in raizes_num.keys() for _ in range(int(raizes_num[z]))]
        self.polos = [complex(p.evalf()) for p in raizes_den.keys() for _ in range(int(raizes_den[p]))]
        self.nz = len(self.zeros)
        self.np = len(self.polos)

        all_reals = [p.real for p in self.polos] + [z.real for z in self.zeros]
        if all_reals:
            self.min_x = min(all_reals)
            self.max_x = max(all_reals)
            self.span = max(self.max_x - self.min_x, 4.0)

        self._add_polos_zeros()

        return {
            "num_expr": self.num_expr,
            "den_expr": self.den_expr,
            "P_expr": self.P_expr,
            "zeros": self.zeros,
            "polos": self.polos,
            "ganho_constante": self.ganho_constante,
        }

    def dados_passo1(self) -> Dict[str, Any]:
        return {
            "gh": self.P_expr,
            "num_expr": self.num_expr,
            "den_expr": self.den_expr,
            "char_expr": self.char_expr,
            "char_derivative": self.char_derivative,
            "char_poly_coeffs": self.char_poly.all_coeffs(),
        }

    # ------------------------------------------------------------------
    # Passo 4: eixo real
    # ------------------------------------------------------------------
    def calcular_segmentos_eixo_real(self) -> Dict[str, Any]:
        pontos_reais = sorted(
            [p.real for p in self.polos if abs(p.imag) < 1e-8]
            + [z.real for z in self.zeros if abs(z.imag) < 1e-8]
        )

        # Mantém multiplicidades. Na contagem de pontos à direita, multiplicidades
        # também devem ser contabilizadas.
        limites = [-math.inf] + pontos_reais + [math.inf]
        segmentos = []
        testes = []

        for esq, dir_ in zip(limites[:-1], limites[1:]):
            if math.isinf(esq) and math.isinf(dir_):
                continue
            if math.isinf(esq):
                teste = dir_ - max(self.span * 0.15, 1.0)
            elif math.isinf(dir_):
                teste = esq + max(self.span * 0.15, 1.0)
            else:
                teste = (esq + dir_) / 2.0

            direita = sum(1 for p in pontos_reais if p > teste)
            pertence = direita % 2 == 1
            testes.append(
                {
                    "esquerda": esq,
                    "direita": dir_,
                    "ponto_teste": teste,
                    "elementos_direita": direita,
                    "paridade": "ímpar" if pertence else "par",
                    "pertence": pertence,
                }
            )
            if pertence:
                segmentos.append((esq, dir_))

        self.real_segments = segmentos

        return {
            "pontos_reais": pontos_reais,
            "testes": testes,
            "segmentos": segmentos,
        }

    # ------------------------------------------------------------------
    # Passos 5, 6 e 7
    # ------------------------------------------------------------------
    def calcular_assintotas(self) -> Dict[str, Any]:
        n_ass = self.np - self.nz
        if n_ass <= 0:
            self.assymptote_center = None
            self.assymptote_angles = []
            return {"numero": n_ass, "centro": None, "angulos": []}

        soma_polos = sum(self.polos)
        soma_zeros = sum(self.zeros)
        centro = (soma_polos - soma_zeros) / n_ass
        centro = float(np.real(centro))
        angulos = [((2 * q + 1) * 180.0) / n_ass for q in range(n_ass)]

        self.assymptote_center = centro
        self.assymptote_angles = angulos

        return {
            "numero": n_ass,
            "soma_polos": soma_polos,
            "soma_zeros": soma_zeros,
            "centro": centro,
            "angulos": angulos,
        }

    # ------------------------------------------------------------------
    # Passo 8: saída/chegada
    # ------------------------------------------------------------------
    def calcular_pontos_saida_chegada(self) -> Dict[str, Any]:
        # K(s) = -D(s)/N(s)
        K_of_s = sp.cancel(-self.den_expr / self.num_expr)
        dK_ds = sp.simplify(sp.diff(K_of_s, self.s))
        numerador_derivada = sp.factor(sp.together(dK_ds).as_numer_denom()[0])
        denominador_derivada = sp.factor(sp.together(dK_ds).as_numer_denom()[1])

        # Forma equivalente sem quociente, útil para a resolução manual.
        Dp = sp.diff(self.den_expr, self.s)
        Np = sp.diff(self.num_expr, self.s)
        equacao_numerador = sp.factor(Dp * self.num_expr - self.den_expr * Np)

        candidatos = []
        try:
            roots = sp.nroots(equacao_numerador)
            for r in roots:
                z = complex(r)
                if abs(z.imag) > 1e-6:
                    K_val = None
                    try:
                        K_val_complex = complex(sp.N(K_of_s.subs(self.s, z)))
                        K_val = K_val_complex.real if abs(K_val_complex.imag) < 1e-6 else K_val_complex
                    except Exception:
                        pass
                    candidatos.append({
                        "s": z,
                        "K": K_val,
                        "real": False,
                        "no_lgr": False,
                        "valido": False,
                    })
                    continue

                s_real = z.real
                K_val_complex = complex(sp.N(K_of_s.subs(self.s, s_real)))
                K_val = K_val_complex.real if abs(K_val_complex.imag) < 1e-6 else K_val_complex
                no_lgr = self._ponto_no_segmento_real(s_real)
                valido = no_lgr and isinstance(K_val, (int, float, np.floating)) and K_val > 0
                candidatos.append({
                    "s": complex(s_real, 0),
                    "K": float(K_val) if isinstance(K_val, (int, float, np.floating)) else K_val,
                    "real": True,
                    "no_lgr": no_lgr,
                    "valido": valido,
                })
        except Exception:
            candidatos = []

        self.breakaway_candidates = candidatos
        return {
            "K_of_s": K_of_s,
            "K_num": self.den_expr,
            "K_den": self.num_expr,
            # Aliases mantidos para compatibilidade com versões anteriores da interface.
            "den_expr": self.den_expr,
            "num_expr": self.num_expr,
            "dK_ds": dK_ds,
            "numerador_derivada": numerador_derivada,
            "denominador_derivada": denominador_derivada,
            "Dp": Dp,
            "Np": Np,
            "equacao_numerador": equacao_numerador,
            "candidatos": candidatos,
        }

    def _ponto_no_segmento_real(self, x: float) -> bool:
        if not self.real_segments:
            self.calcular_segmentos_eixo_real()
        for a, b in self.real_segments:
            left_ok = True if math.isinf(a) else x > a + 1e-7
            right_ok = True if math.isinf(b) else x < b - 1e-7
            if left_ok and right_ok:
                return True
        return False

    # ------------------------------------------------------------------
    # Passo 9: Routh-Hurwitz
    # ------------------------------------------------------------------
    def routh_table(self) -> List[Tuple[sp.Expr, List[sp.Expr]]]:
        coeffs = list(self.char_poly.all_coeffs())
        n = len(coeffs) - 1
        cols = int(math.ceil((n + 1) / 2))
        table = [[sp.Integer(0) for _ in range(cols)] for _ in range(n + 1)]

        row_powers = list(range(n, -1, -1))
        rows = [s for s in row_powers]

        table[0][: len(coeffs[0::2])] = coeffs[0::2]
        table[1][: len(coeffs[1::2])] = coeffs[1::2]

        epsilon = sp.Symbol(r"\epsilon", positive=True)
        used_epsilon = False
        for i in range(2, n + 1):
            for j in range(cols - 1):
                a = table[i - 2][0]
                b = table[i - 2][j + 1]
                c = table[i - 1][0]
                d = table[i - 1][j + 1]
                if sp.simplify(c) == 0:
                    c = epsilon
                    used_epsilon = True
                table[i][j] = sp.factor((c * b - a * d) / c)

        # O indicador é retornado em atributo auxiliar para a interface poder
        # informar que o caso especial epsilon foi necessário.
        self.routh_used_epsilon = used_epsilon

        return [(sp.Integer(rows[i]), [sp.factor(v) for v in table[i]]) for i in range(n + 1)]

    def calcular_routh(self) -> Dict[str, Any]:
        table = self.routh_table()
        first_column = [row[1][0] for row in table]

        candidates = []
        for expr in first_column[1:]:
            if expr in (0, sp.nan):
                continue
            try:
                sols = sp.solve(sp.Eq(sp.together(expr), 0), self.K)
            except Exception:
                sols = []
            for sol in sols:
                sol_num = complex(sp.N(sol))
                if abs(sol_num.imag) < 1e-8 and sol_num.real > 0:
                    candidates.append(float(sol_num.real))

        # Remove duplicatas numéricas.
        K_candidates = []
        for value in sorted(candidates):
            if not any(abs(value - old) < 1e-7 for old in K_candidates):
                K_candidates.append(value)

        detailed = []
        for kval in K_candidates:
            char_at_k = sp.Poly(sp.expand(self.char_expr.subs(self.K, kval)), self.s)
            roots = np.roots([float(c) for c in char_at_k.all_coeffs()])
            imag_roots = [r for r in roots if abs(r.real) < 1e-4]

            aux_info = self._auxiliary_polynomial_from_routh(kval, table)
            detailed.append({
                "K": kval,
                "roots": roots,
                "imaginary_roots": imag_roots,
                "auxiliary": aux_info,
            })

        return {
            "table": table,
            "first_column": first_column,
            "K_candidates": K_candidates,
            "crossing_candidates": detailed,
            "used_epsilon": self.routh_used_epsilon,
        }

    def _auxiliary_polynomial_from_routh(self, K_value: float, table):
        # Localiza a linha que zera integralmente após a substituição de K.
        evaluated = []
        for power, row in table:
            vals = [sp.N(expr.subs(self.K, K_value)) if expr not in (sp.nan,) else sp.nan for expr in row]
            evaluated.append((int(power), vals))

        zero_row_index = None
        for i, (_power, vals) in enumerate(evaluated):
            finite_vals = [v for v in vals if v != sp.nan]
            if finite_vals and all(abs(complex(v)) < 1e-7 for v in finite_vals):
                if i >= 2:
                    zero_row_index = i
                    break

        if zero_row_index is None:
            # Também podemos construir o polinômio auxiliar a partir das raízes
            # puramente imaginárias, mas sem a estrutura de Routh explícita.
            return None

        source_power, source_vals = evaluated[zero_row_index - 1]
        terms = []
        for j, value in enumerate(source_vals):
            if value == sp.nan:
                continue
            exponent = source_power - 2 * j
            if exponent < 0:
                break
            if abs(complex(value)) > 1e-12:
                terms.append((exponent, float(sp.N(value))))
        poly = sum(c * self.s ** e for e, c in terms)
        if terms:
            max_exp = max(e for e, _c in terms)
            coeffs = [0.0] * (max_exp + 1)
            for exponent, coeff in terms:
                coeffs[max_exp - exponent] = coeff
            aux_roots = np.roots(coeffs)
        else:
            aux_roots = []
        return {
            "row_power": source_power,
            "terms": terms,
            "polynomial": sp.expand(poly),
            "roots": aux_roots,
        }

    def _cruzamento_por_equacao_s2(self, routh: Dict[str, Any]) -> Dict[str, Any]:
        """Aplica explicitamente a equação da linha s² para obter K no cruzamento.

        Para o caso de quarta ordem do material, a linha s² tem a forma
        b1*s² + b2 = 0. Com s=j*w: -b1*w² + b2 = 0.
        O valor de w é obtido pela parte imaginária da equação característica e,
        em seguida, K é obtido substituindo esse w na equação auxiliar da linha s².
        A equação real da característica é usada como verificação final.
        """
        result = {
            "disponivel": False,
            "linha_s2": None,
            "auxiliary": None,
            "aux_jw": None,
            "char_re": None,
            "char_im": None,
            "omega_equation": None,
            "k_expression": None,
            "omega_values": [],
            "calculos": [],
            "cruzamentos": [],
        }

        row_s2 = None
        for power, row in routh["table"]:
            if int(power) == 2:
                row_s2 = row
                break
        if row_s2 is None or len(row_s2) < 2:
            return result

        a = sp.factor(sp.sympify(row_s2[0]))
        b = sp.factor(sp.sympify(row_s2[1]))
        if a == 0 and b == 0:
            return result

        aux = sp.expand(a * self.s**2 + b)
        w = self.w
        aux_jw = sp.factor(sp.expand(aux.subs(self.s, sp.I * w)).as_real_imag()[0])
        try:
            k_solutions = sp.solve(sp.Eq(aux_jw, 0), self.K)
            k_expression = sp.factor(k_solutions[0]) if k_solutions else None
        except Exception:
            k_expression = None
        # Equação característica em s=j*w.
        phi_jw = sp.expand(self.char_expr.subs(self.s, sp.I * w))
        char_re = sp.factor(sp.re(phi_jw).expand())
        char_im = sp.factor(sp.im(phi_jw).expand())

        # Para w != 0, a parte imaginária é dividida por w. Isso reproduz a
        # substituição usada manualmente no caso típico de ordem 4.
        if sp.simplify(char_im.subs(w, 0)) == 0:
            omega_equation = sp.factor(sp.cancel(char_im / w))
        else:
            omega_equation = char_im

        omega_values = []
        try:
            roots_omega = sp.nroots(sp.Poly(omega_equation, w))
            for root in roots_omega:
                z = complex(root)
                if abs(z.imag) < 1e-7 and z.real > 1e-7:
                    omega_values.append(float(z.real))
        except Exception:
            pass

        # Fallback simbólico para equações simples em w².
        if not omega_values:
            try:
                sols = sp.solve(sp.Eq(omega_equation, 0), w)
                for sol in sols:
                    z = complex(sp.N(sol))
                    if abs(z.imag) < 1e-7 and z.real > 1e-7:
                        omega_values.append(float(z.real))
            except Exception:
                pass

        calculations = []
        crossings = []
        for omega in omega_values:
            aux_num = sp.N(aux_jw.subs(w, omega))
            try:
                K_solutions = sp.solve(sp.Eq(aux_num, 0), self.K)
            except Exception:
                K_solutions = []
            for K_sol in K_solutions:
                try:
                    kval = complex(sp.N(K_sol))
                except Exception:
                    continue
                if abs(kval.imag) > 1e-7 or kval.real <= 0:
                    continue
                kval_f = float(kval.real)
                real_residual = float(sp.N(char_re.subs({w: omega, self.K: kval_f})))
                imag_residual = float(sp.N(char_im.subs({w: omega, self.K: kval_f})))
                valid = abs(real_residual) < 1e-4 and abs(imag_residual) < 1e-4
                calculations.append({
                    "w": omega,
                    "K": kval_f,
                    "s": complex(0.0, omega),
                    "real_residual": real_residual,
                    "imag_residual": imag_residual,
                    "valido": valid,
                })
                if valid:
                    crossings.append({"w": omega, "K": kval_f, "origem": "Equação de s²"})

        # Evita duplicatas.
        unique = []
        for c in crossings:
            if not any(abs(c["w"] - u["w"]) < 1e-5 and abs(c["K"] - u["K"]) < 1e-5 for u in unique):
                unique.append(c)

        result.update({
            "disponivel": True,
            "linha_s2": [a, b],
            "auxiliary": aux,
            "aux_jw": aux_jw,
            "char_re": char_re,
            "char_im": char_im,
            "omega_equation": omega_equation,
            "k_expression": k_expression,
            "omega_values": omega_values,
            "calculos": calculations,
            "cruzamentos": unique,
        })
        return result

    # ------------------------------------------------------------------
    # Passo 9 auxiliar: cruzamento analítico
    # ------------------------------------------------------------------
    def calcular_cruzamento_eixo_imaginario(self) -> Dict[str, Any]:
        # Método geral usando s = j*w. Re e Im da equação D(jw) + K N(jw)=0.
        D_jw = sp.expand(self.den_expr.subs(self.s, sp.I * self.w))
        N_jw = sp.expand(self.num_expr.subs(self.s, sp.I * self.w))

        re_D, im_D = sp.re(D_jw).expand(), sp.im(D_jw).expand()
        re_N, im_N = sp.re(N_jw).expand(), sp.im(N_jw).expand()

        # Eliminando K entre as partes real e imaginária.
        eliminacao = sp.factor(im_D * re_N - re_D * im_N)
        candidatos_w = []
        cruzamentos = []

        try:
            pol_w = sp.Poly(eliminacao, self.w)
            for root in sp.nroots(pol_w):
                root_c = complex(root)
                if abs(root_c.imag) > 1e-7:
                    continue
                w_val = float(root_c.real)
                if w_val <= 1e-7:
                    continue

                # Calcula K a partir de qualquer componente não degenerada.
                k_val = None
                re_n = float(sp.N(re_N.subs(self.w, w_val)))
                im_n = float(sp.N(im_N.subs(self.w, w_val)))
                re_d = float(sp.N(re_D.subs(self.w, w_val)))
                im_d = float(sp.N(im_D.subs(self.w, w_val)))

                if abs(re_n) > 1e-8:
                    k_val = -re_d / re_n
                elif abs(im_n) > 1e-8:
                    k_val = -im_d / im_n

                candidatos_w.append(w_val)
                if k_val is not None and k_val > 0:
                    cruzamentos.append({"w": w_val, "K": float(k_val)})
        except Exception:
            pass

        # Método solicitado para a resolução manual: usar a equação da linha s².
        routh = self.calcular_routh()
        s2_method = self._cruzamento_por_equacao_s2(routh)
        for item in s2_method["cruzamentos"]:
            if not any(abs(item["w"] - c["w"]) < 1e-4 and abs(item["K"] - c["K"]) < 1e-4 for c in cruzamentos):
                cruzamentos.append(item)

        # União com os pontos obtidos diretamente da análise das raízes do Routh.
        for item in routh["crossing_candidates"]:
            kval = item["K"]
            for rr in item["imaginary_roots"]:
                if rr.imag > 1e-6:
                    candidato = {"w": float(rr.imag), "K": float(kval), "origem": "Routh-Hurwitz"}
                    if not any(abs(candidato["w"] - c["w"]) < 1e-4 and abs(candidato["K"] - c["K"]) < 1e-4 for c in cruzamentos):
                        cruzamentos.append(candidato)

        self.crossings = cruzamentos
        return {
            "D_jw": D_jw,
            "N_jw": N_jw,
            "re_D": re_D,
            "im_D": im_D,
            "re_N": re_N,
            "im_N": im_N,
            "eliminacao": eliminacao,
            "candidatos_w": candidatos_w,
            "cruzamentos": cruzamentos,
            "routh": routh,
            "s2_method": s2_method,
        }

    # ------------------------------------------------------------------
    # Passo 10: ângulos de partida/chegada
    # ------------------------------------------------------------------
    def analisar_angulo_partida_chegada(self) -> Dict[str, Any]:
        resultados = []

        for p in self.polos:
            if abs(p.imag) < 1e-7:
                continue
            vet_polos = self._vetores_a_partir_de(p, [q for q in self.polos if abs(q - p) > 1e-7])
            vet_zeros = self._vetores_a_partir_de(p, self.zeros)
            soma_theta = sum(v.angulo_deg for v in vet_polos)
            soma_phi = sum(v.angulo_deg for v in vet_zeros)
            ang = self._normalize_angle_360(180.0 - soma_theta + soma_phi)
            resultados.append({
                "tipo": "partida",
                "ponto": p,
                "vetores_polos": vet_polos,
                "vetores_zeros": vet_zeros,
                "soma_polos": soma_theta,
                "soma_zeros": soma_phi,
                "angulo": ang,
            })

        for z in self.zeros:
            if abs(z.imag) < 1e-7:
                continue
            vet_zeros = self._vetores_a_partir_de(z, [q for q in self.zeros if abs(q - z) > 1e-7])
            vet_polos = self._vetores_a_partir_de(z, self.polos)
            soma_phi = sum(v.angulo_deg for v in vet_zeros)
            soma_theta = sum(v.angulo_deg for v in vet_polos)
            ang = self._normalize_angle_360(180.0 - soma_phi + soma_theta)
            resultados.append({
                "tipo": "chegada",
                "ponto": z,
                "vetores_polos": vet_polos,
                "vetores_zeros": vet_zeros,
                "soma_polos": soma_theta,
                "soma_zeros": soma_phi,
                "angulo": ang,
            })

        self.departure_arrival = resultados
        return {"resultados": resultados}

    def _vetores_a_partir_de(self, destino: complex, origens: List[complex]) -> List[VetorAnalise]:
        result = []
        for origem in origens:
            vetor = destino - origem
            result.append(
                VetorAnalise(
                    origem=origem,
                    destino=destino,
                    dx=float(vetor.real),
                    dy=float(vetor.imag),
                    magnitude=float(abs(vetor)),
                    angulo_deg=self._angle_deg(vetor),
                )
            )
        return result

    # ------------------------------------------------------------------
    # Passos 11 e 12: ponto de teste
    # ------------------------------------------------------------------
    def analisar_ponto_teste(self, s_teste: complex) -> Dict[str, Any]:
        s_teste = complex(s_teste)
        vet_polos = self._vetores_a_partir_de(s_teste, self.polos)
        vet_zeros = self._vetores_a_partir_de(s_teste, self.zeros)

        soma_polos = sum(v.angulo_deg for v in vet_polos)
        soma_zeros = sum(v.angulo_deg for v in vet_zeros)

        # P(s) = C * prod(s-z_i) / prod(s-p_i).
        # A fase de C é 0° para C>0 e 180° para C<0.
        fase_bruta = self.fase_constante_deg + soma_zeros - soma_polos
        fase_mod = fase_bruta % 360.0
        erro_angulo = abs(((fase_mod - 180.0 + 180.0) % 360.0) - 180.0)
        pertence = erro_angulo < 1e-2

        mod_polos = float(np.prod([v.magnitude for v in vet_polos])) if vet_polos else 1.0
        mod_zeros = float(np.prod([v.magnitude for v in vet_zeros])) if vet_zeros else 1.0
        ganho_abs = abs(float(sp.N(self.ganho_constante)))
        # |K C| prod|s-z|/prod|s-p| = 1.
        K_val = mod_polos / (ganho_abs * mod_zeros) if pertence else None

        return {
            "s_teste": s_teste,
            "vetores_polos": vet_polos,
            "vetores_zeros": vet_zeros,
            "soma_angulos_polos": soma_polos,
            "soma_angulos_zeros": soma_zeros,
            "fase_constante_deg": self.fase_constante_deg,
            "ganho_constante": self.ganho_constante,
            "fase_bruta": fase_bruta,
            "fase_mod": fase_mod,
            "erro_angulo": erro_angulo,
            "pertence": pertence,
            "produto_mod_polos": mod_polos,
            "produto_mod_zeros": mod_zeros,
            "K": K_val,
        }

    # ------------------------------------------------------------------
    # Gráfico exato
    # ------------------------------------------------------------------
    def adicionar_elementos_geometricos(self):
        # Buffer invisível apenas para preservar a escala gráfica.
        self.fig.add_trace(
            go.Scatter(
                x=[self.min_x - self.span * 0.10, self.max_x + self.span * 0.10],
                y=[-self.span * 0.5, self.span * 0.5],
                mode="markers",
                marker=dict(color="rgba(0,0,0,0)"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

        raio = self.span * 1.5
        if self.assymptote_center is not None:
            for i, ang in enumerate(self.assymptote_angles):
                rad = math.radians(ang)
                dx = 0.0 if abs(math.cos(rad)) < 1e-10 else math.cos(rad)
                dy = 0.0 if abs(math.sin(rad)) < 1e-10 else math.sin(rad)
                self.fig.add_trace(
                    go.Scatter(
                        x=[self.assymptote_center, self.assymptote_center + raio * dx],
                        y=[0, raio * dy],
                        mode="lines",
                        line=dict(color="gray", width=1.5, dash="dash"),
                        name=f"Assíntota {i + 1}",
                    )
                )

        for i, (a, b) in enumerate(self.real_segments):
            x_start = a if not math.isinf(a) else b - self.span
            x_end = b if not math.isinf(b) else a + self.span
            self.fig.add_trace(
                go.Scatter(
                    x=[x_start, x_end],
                    y=[0, 0],
                    mode="lines",
                    line=dict(color="green", width=4),
                    name=f"Segmento Real {i + 1}",
                )
            )

    def adicionar_marcadores_especiais(self):
        validos = [c for c in self.breakaway_candidates if c.get("valido") and c.get("real")]
        if validos:
            xs = [c["s"].real for c in validos]
            self.fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=[0] * len(xs),
                    mode="markers",
                    marker=dict(symbol="square", size=9, color="purple"),
                    name="Saída/Chegada",
                )
            )

        if self.crossings:
            xs = [0 for _ in self.crossings for _ in (0, 1)]
            ys = [v for c in self.crossings for v in (c["w"], -c["w"])]
            self.fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="markers",
                    marker=dict(symbol="diamond", size=9, color="orange"),
                    name="Cruzamento Im",
                )
            )

    def adicionar_ponto_teste(self, s_teste: complex):
        self.fig.add_trace(
            go.Scatter(
                x=[s_teste.real],
                y=[s_teste.imag],
                mode="markers",
                marker=dict(symbol="star", size=14, color="magenta", line=dict(width=1, color="black")),
                name="Ponto Teste",
            )
        )

    def _calcular_lgr_numerico(self):
        """Calcula numericamente os ramos do LGR sem depender do python-control.

        Resolve diretamente o polinômio característico

            D(s) + K N(s) = 0

        para uma sequência de valores de K. As raízes são então ordenadas entre
        pontos consecutivos pelo menor custo de deslocamento. O método é escrito
        para ser robusto a polos múltiplos e a eventuais valores numéricos
        inválidos retornados por ``numpy.roots``.
        """
        # K=0 é incluído explicitamente para que os polos de malha aberta apareçam
        # no início do LGR. O restante usa escala logarítmica para acompanhar a
        # evolução dos ramos até os zeros finitos/assíntotas.
        klist = np.concatenate(([0.0], np.logspace(-5, 5, 1000)))
        roots_all = []
        expected_degree = self.char_poly.degree()

        for kval in klist:
            try:
                poly = sp.Poly(
                    sp.expand(self.char_expr.subs(self.K, float(kval))),
                    self.s,
                )
                coeffs = np.asarray(
                    [complex(c.evalf()) for c in poly.all_coeffs()],
                    dtype=complex,
                )

                # Remove apenas coeficientes líderes efetivamente nulos. Isso
                # evita mudar artificialmente o grau por ruído numérico.
                while len(coeffs) > 1 and abs(coeffs[0]) < 1e-14:
                    coeffs = coeffs[1:]

                if len(coeffs) < 2 or not np.isfinite(coeffs).all():
                    roots = np.full(expected_degree, np.nan + 1j * np.nan, dtype=complex)
                else:
                    roots = np.asarray(np.roots(coeffs), dtype=complex)
                    if len(roots) != expected_degree or not np.isfinite(roots).all():
                        roots = np.full(expected_degree, np.nan + 1j * np.nan, dtype=complex)
            except Exception:
                roots = np.full(expected_degree, np.nan + 1j * np.nan, dtype=complex)

            roots_all.append(roots)

        # O grau esperado é constante para K finito. Mantemos uma matriz
        # retangular e marcamos amostras inválidas com NaN; elas são tratadas
        # explicitamente no acompanhamento dos ramos e não são enviadas ao
        # algoritmo de atribuição.
        nbranches = expected_degree
        ordered = np.full((len(klist), nbranches), np.nan + 1j * np.nan, dtype=complex)

        # Primeira amostra: os polos de malha aberta. Se ela for válida, usamos
        # diretamente; caso contrário, ordenamos deterministicamente.
        first = roots_all[0]
        if np.isfinite(first).all() and len(first) == nbranches:
            ordered[0] = first
        else:
            valid = [r for r in first if np.isfinite(r)]
            ordered[0, :len(valid)] = np.asarray(valid, dtype=complex)

        def assign_roots(previous, current):
            """Associa raízes válidas sem permitir NaN/inf no custo."""
            prev_valid = [i for i, z in enumerate(previous) if np.isfinite(z)]
            curr_valid = [j for j, z in enumerate(current) if np.isfinite(z)]

            result = np.full_like(current, np.nan + 1j * np.nan)
            if not prev_valid:
                if curr_valid:
                    # Ordem determinística para uma primeira amostra recuperada.
                    for pos, j in enumerate(curr_valid):
                        result[pos] = current[j]
                return result

            if not curr_valid:
                return result

            # Há poucos ramos no tipo de exercício da disciplina. Para n <= 8,
            # testar permutações fornece a mesma ideia do problema de atribuição
            # sem depender de scipy e sem aceitar custos inválidos.
            if len(prev_valid) <= 8 and len(curr_valid) == len(prev_valid):
                import itertools

                best_perm = None
                best_cost = float("inf")
                prev_values = [previous[i] for i in prev_valid]
                curr_values = [current[j] for j in curr_valid]
                for perm in itertools.permutations(range(len(curr_values))):
                    cost = sum(abs(prev_values[a] - curr_values[perm[a]]) for a in range(len(prev_values)))
                    if np.isfinite(cost) and cost < best_cost:
                        best_cost = float(cost)
                        best_perm = perm

                if best_perm is not None:
                    for a, prev_index in enumerate(prev_valid):
                        result[prev_index] = curr_values[best_perm[a]]
                    return result

            # Fallback guloso para ordens maiores ou conjuntos com cardinalidade
            # diferente. Nunca calcula distância com entrada inválida.
            unused = set(curr_valid)
            for prev_index in prev_valid:
                if not unused:
                    break
                prev_value = previous[prev_index]
                best_index = min(
                    unused,
                    key=lambda idx: abs(prev_value - current[idx]),
                )
                if np.isfinite(prev_value) and np.isfinite(current[best_index]):
                    result[prev_index] = current[best_index]
                unused.remove(best_index)

            # Se sobraram raízes atuais (situação patológica), elas são colocadas
            # nas posições ainda vazias, sem tentar criar correspondências inválidas.
            free_positions = [i for i, z in enumerate(result) if not np.isfinite(z)]
            for pos, idx in zip(free_positions, sorted(unused)):
                result[pos] = current[idx]
            return result

        for i in range(1, len(klist)):
            ordered[i] = assign_roots(ordered[i - 1], roots_all[i])

        return ordered, klist

    def calcular_lgr_exato(self):
        """Traça numericamente os ramos do LGR.

        Tenta utilizar ``python-control`` por compatibilidade com as versões
        anteriores. Entretanto, algumas versões dessa biblioteca podem falhar
        em ``root_locus`` quando há polos múltiplos ou raízes que se aproximam
        muito entre si, produzindo uma exceção de dimensões no ``vstack``.
        Nesse caso, o método faz automaticamente o cálculo direto de
        D(s) + K N(s) = 0, sem interromper a resolução.
        """
        rlist = None
        klist = None

        if ctrl is not None:
            try:
                rlist, klist = ctrl.root_locus(self.GH, plot=False)

                # Verifica se a biblioteca realmente devolveu uma matriz
                # retangular válida para o gráfico.
                rlist = np.asarray(rlist, dtype=complex)
                klist = np.asarray(klist, dtype=float)
                if rlist.ndim != 2 or klist.ndim != 1 or rlist.shape[0] != len(klist):
                    raise ValueError("root_locus retornou dimensões incompatíveis.")
            except Exception:
                rlist, klist = self._calcular_lgr_numerico()
        else:
            rlist, klist = self._calcular_lgr_numerico()

        self._last_rlist = rlist
        self._last_klist = klist

        for i in range(rlist.shape[1]):
            self.fig.add_trace(
                go.Scatter(
                    x=np.real(rlist[:, i]),
                    y=np.imag(rlist[:, i]),
                    mode="lines",
                    line=dict(width=2.5),
                    name=f"Ramo Final {i + 1}",
                )
            )
        return rlist, klist
