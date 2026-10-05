import math
import os
import shutil
import subprocess
import tempfile
import warnings
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

try:
    import control as ctrl
except ImportError:  # Permite executar os cálculos mesmo sem python-control instalado.
    ctrl = None
import numpy as np
import sympy as sp
import plotly.graph_objects as go

from io import BytesIO
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image,
)

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

    # ------------------------------------------------------------------
    # Exportação da resolução em PDF
    # ------------------------------------------------------------------
    @staticmethod
    def _pdf_find_font(weight="normal"):
        """Localiza uma fonte TTF de forma portável em Windows, Linux e macOS.

        Prioriza as fontes distribuídas com o Matplotlib, que normalmente estão
        disponíveis junto da própria instalação do Python e não dependem de
        caminhos específicos do sistema operacional.
        """
        candidates = []

        # 1) Fontes do próprio Matplotlib (preferência).
        try:
            prop = fm.FontProperties(family="DejaVu Sans", weight=weight)
            found = fm.findfont(prop, fallback_to_default=True)
            if found:
                candidates.append(Path(found))
        except Exception:
            pass

        # 2) Caminhos comuns do sistema como fallback.
        if weight == "bold":
            names = [
                "DejaVuSans-Bold.ttf",
                "NotoSans-Bold.ttf",
                "LiberationSans-Bold.ttf",
                "arialbd.ttf",
                "segoeuib.ttf",
            ]
        else:
            names = [
                "DejaVuSans.ttf",
                "NotoSans-Regular.ttf",
                "LiberationSans-Regular.ttf",
                "arial.ttf",
                "segoeui.ttf",
            ]

        system_dirs = [
            Path("/usr/share/fonts/truetype/dejavu"),
            Path("/usr/share/fonts/truetype/noto"),
            Path("/usr/share/fonts/truetype/liberation"),
            Path("/usr/local/share/fonts"),
            Path("C:/Windows/Fonts"),
            Path("C:/Windows/Fonts/Arial"),
            Path("/System/Library/Fonts"),
            Path("/Library/Fonts"),
        ]

        for directory in system_dirs:
            for name in names:
                candidates.append(directory / name)

        for candidate in candidates:
            try:
                if candidate.is_file():
                    return str(candidate)
            except OSError:
                continue

        return None

    @staticmethod
    def _pdf_font_setup():
        """Configura fontes PDF sem depender de um caminho Linux fixo.

        Em instalações normais, usa DejaVu Sans fornecida pelo Matplotlib.
        Caso nenhuma TTF seja encontrada, usa as fontes padrão do ReportLab
        (Helvetica/Helvetica-Bold), evitando que a exportação quebre por causa
        de uma fonte ausente.
        """
        regular = AnalisadorLGR._pdf_find_font("normal")
        bold = AnalisadorLGR._pdf_find_font("bold")

        if regular and bold:
            if "LGR-Regular" not in pdfmetrics.getRegisteredFontNames():
                pdfmetrics.registerFont(TTFont("LGR-Regular", regular))
            if "LGR-Bold" not in pdfmetrics.getRegisteredFontNames():
                pdfmetrics.registerFont(TTFont("LGR-Bold", bold))
            return "LGR-Regular", "LGR-Bold"

        # Último fallback: fontes base 14 do PDF, que não exigem arquivos TTF.
        return "Helvetica", "Helvetica-Bold"

    @staticmethod
    def _pdf_expr(expr):
        """Representação legível de expressões SymPy em PDF textual."""
        if expr is None:
            return "-"
        try:
            text = str(sp.expand(expr))
        except Exception:
            text = str(expr)
        text = text.replace("**", "^")
        text = text.replace("*", "·")
        text = text.replace("sqrt", "√")
        text = text.replace("I", "j")
        return text

    @staticmethod
    def _tex_expr(expr):
        """Renderiza uma expressão SymPy como LaTeX, tentando eliminar ruído numérico."""
        if expr is None:
            return r"-"
        try:
            clean = sp.nsimplify(sp.expand(expr), tolerance=1e-8, rational=True)
            clean = sp.cancel(clean)
            return sp.latex(clean)
        except Exception:
            return sp.latex(sp.expand(sp.sympify(expr)))

    @staticmethod
    def _tex_escape(text):
        """Escapa texto comum para uso dentro do documento LaTeX."""
        replacements = {
            "\\": r"\textbackslash{}",
            "&": r"\&",
            "%": r"\%",
            "$": r"\$",
            "#": r"\#",
            "_": r"\_",
            "{": r"\{",
            "}": r"\}",
            "~": r"\textasciitilde{}",
            "^": r"\textasciicircum{}",
        }
        return "".join(replacements.get(ch, ch) for ch in str(text))

    @staticmethod
    def _tex_complex(z, digits=4):
        """Forma compacta de número complexo para matemática LaTeX."""
        z = complex(z)
        if abs(z.imag) < 10 ** (-digits):
            return f"{z.real:.{digits}f}"
        sign = "+" if z.imag >= 0 else "-"
        return rf"{z.real:.{digits}f} {sign} {abs(z.imag):.{digits}f}j"

    @staticmethod
    def _tex_interval(a, b, digits=4):
        left = r"-\infty" if math.isinf(a) and a < 0 else f"{a:.{digits}f}"
        right = r"+\infty" if math.isinf(b) and b > 0 else f"{b:.{digits}f}"
        return rf"({left}, {right})"

    def gerar_pdf_modo_prova(
        self, p1, p4, p7, p8, p9, p10, ptest, ponto_teste,
    ) -> bytes:
        """Gera uma versão enxuta do PDF, otimizada para copiar durante a prova.

        O documento é compilado com LaTeX real (pdflatex), para melhorar a
        leitura de fórmulas, frações, subscritos e símbolos matemáticos.
        """
        self._pdf_test_point = ponto_teste
        if self._last_rlist is None:
            self.calcular_lgr_exato()

        # O LaTeX é o caminho preferencial. Se não estiver instalado, usamos
        # um fallback COMPACTO próprio do Modo Prova -- nunca o PDF completo.
        pdflatex = shutil.which("pdflatex")
        if pdflatex is None:
            return self._gerar_pdf_modo_prova_fallback(
                p1, p4, p7, p8, p9, p10, ptest, ponto_teste
            )

        def eq(expr):
            return rf"\[{expr}\]"

        doc_lines = [
            r"\documentclass[10pt,a4paper]{article}",
            r"\usepackage[utf8]{inputenc}",
            r"\usepackage[T1]{fontenc}",
            r"\usepackage[brazil]{babel}",
            r"\usepackage{lmodern}",
            r"\usepackage{amsmath,amssymb,mathtools}",
            r"\usepackage[a4paper,margin=1.35cm]{geometry}",
            r"\usepackage{graphicx}",
            r"\usepackage{booktabs}",
            r"\usepackage{microtype}",
            r"\setlength{\parindent}{0pt}",
            r"\setlength{\parskip}{3pt}",
            r"\newcommand{\passo}[1]{\vspace{7pt}\par\noindent\textbf{#1}\par\smallskip\nobreak}",
            r"\newcommand{\resultado}[1]{\vspace{1pt}\fbox{\parbox{0.97\linewidth}{#1}}\vspace{2pt}}",
            r"\begin{document}",
            r"\begin{center}",
            r"{\Large\bfseries Lugar Geométrico das Raízes -- Modo Prova}\\[-1mm]",
            r"{\small Roteiro enxuto para reprodução manual em papel}",
            r"\end{center}",
        ]

        # Dados do problema
        doc_lines.append(r"\passo{Dados do problema}")
        doc_lines.append(eq(rf"G(s)=\frac{{{self._tex_expr(self._poly_to_expr(self.numG))}}}{{{self._tex_expr(self._poly_to_expr(self.denG))}}}\qquad H(s)=\frac{{{self._tex_expr(self._poly_to_expr(self.numH))}}}{{{self._tex_expr(self._poly_to_expr(self.denH))}}}"))
        doc_lines.append(eq(rf"s_t={self._tex_complex(ponto_teste,4)}"))

        # Passo 1
        doc_lines.append(r"\passo{1. Polinômio característico}")
        doc_lines.append(eq(r"1+K P(s)=0"))
        doc_lines.append(eq(rf"\Phi(s,K)={self._tex_expr(p1['char_expr'])}=0"))

        # Passo 2
        doc_lines.append(r"\passo{2. Pólos e zeros}")
        poles = r"\; ;\;".join(self._tex_complex(x,4) for x in self.polos) if self.polos else r"\varnothing"
        zeros = r"\; ;\;".join(self._tex_complex(x,4) for x in self.zeros) if self.zeros else r"\varnothing"
        doc_lines.append(eq(rf"p_i=\{{{poles}\}},\qquad z_j=\{{{zeros}\}}"))
        doc_lines.append(eq(rf"n_P={self.np},\qquad n_Z={self.nz}"))

        # Passos 3-6
        doc_lines.append(r"\passo{3--6. Esboço, eixo real, ramos e simetria}")
        segs = p4.get("segmentos", [])
        seg_text = r"\; ;\;".join(self._tex_interval(a,b,4) for a,b in segs) if segs else r"\varnothing"
        doc_lines.append(eq(rf"\text{{Segmentos reais do LGR: }}{seg_text}"))
        doc_lines.append(eq(rf"n_{{ramos}}=n_P={self.np}"))
        doc_lines.append(r"\textit{Como os coeficientes são reais, o LGR é simétrico em relação ao eixo real.}")

        # Passo 7
        doc_lines.append(r"\passo{7. Assíntotas}")
        if p7.get("numero", 0) > 0:
            angles = r"\; ;\;".join(f"{a:.2f}^\\circ" for a in p7.get("angulos", []))
            doc_lines.append(eq(rf"q=n_P-n_Z={p7['numero']},\qquad \sigma_A={p7['centro']:.6f}"))
            doc_lines.append(eq(rf"\phi_A={angles}"))
        else:
            doc_lines.append(eq(r"q=n_P-n_Z\le 0\Rightarrow\text{ sem assíntotas para }\infty"))

        # Passo 8
        doc_lines.append(r"\passo{8. Pontos de saída/chegada}")
        doc_lines.append(eq(r"K(s)=-\frac{D(s)}{N(s)},\qquad D'(s)N(s)-D(s)N'(s)=0"))
        valid = [c for c in p8.get("candidatos", []) if c.get("real") and c.get("valido")]
        if valid:
            for c in valid:
                k = c.get("K")
                kval = "--" if k is None else f"{k:.6f}"
                doc_lines.append(eq(rf"s_b={c['s'].real:.6f},\qquad K_b={kval}"))
        else:
            doc_lines.append(r"Não foi encontrado ponto válido para $K>0$ nos segmentos reais.")

        # Passo 9
        doc_lines.append(r"\passo{9. Cruzamento do eixo imaginário (Routh-Hurwitz)}")
        routh = p9.get("routh", {})
        table = routh.get("table", [])
        if table:
            max_cols = max(len(row) for _, row in table)
            cols = "c" * (max_cols + 1)
            doc_lines.append(rf"\begin{{center}}\small\begin{{tabular}}{{{cols}}}")
            doc_lines.append(r"\toprule")
            header = " & ".join([r"Linha"] + [rf"$c_{i+1}$" for i in range(max_cols)]) + " \\\\" 
            doc_lines.append(header)
            doc_lines.append(r"\midrule")
            for power, row in table:
                cells = [rf"$s^{{{int(power)}}}$"] + [rf"${self._tex_expr(v)}$" for v in row]
                cells += [r"$0$"] * (max_cols - len(row))
                doc_lines.append(" & ".join(cells) + " \\\\" )
            doc_lines.append(r"\bottomrule")
            doc_lines.append(r"\end{tabular}\end{center}")
        crossings = p9.get("cruzamentos", [])
        if crossings:
            for c in crossings:
                doc_lines.append(eq(rf"\boxed{{s=\pm j{c['w']:.6f},\qquad K={c['K']:.6f}}}"))
        else:
            doc_lines.append(eq(r"\text{Não há cruzamento detectado para }K>0."))

        # Passo 10
        doc_lines.append(r"\passo{10. Ângulos de partida/chegada}")
        if p10.get("resultados"):
            for item in p10["resultados"]:
                ponto = self._tex_complex(item["ponto"], 4)
                if item["tipo"] == "partida":
                    expr = rf"\theta_{{partida}}=180^\circ-({item['soma_polos']:.4f}^\circ)+({item['soma_zeros']:.4f}^\circ)={item['angulo']:.4f}^\circ"
                else:
                    expr = rf"\theta_{{chegada}}=180^\circ-({item['soma_zeros']:.4f}^\circ)+({item['soma_polos']:.4f}^\circ)={item['angulo']:.4f}^\circ"
                doc_lines.append(eq(rf"s={ponto}:\quad {expr}"))
        else:
            doc_lines.append(r"Não há pólos/zeros complexos que exijam cálculo de ângulo.")

        # Passos 11-12
        doc_lines.append(r"\passo{11. Teste da condição de ângulo}")
        doc_lines.append(eq(rf"\angle P(s_t)={ptest['fase_mod']:.4f}^\circ\qquad(\text{{alvo: }}180^\circ\;\bmod\;360^\circ)"))
        pertence = bool(ptest.get("pertence"))
        status = "PERTENCE" if pertence else "NÃO PERTENCE"
        doc_lines.append(rf"\resultado{{\textbf{{{status} ao LGR para }}$K>0$.}}")

        doc_lines.append(r"\passo{12. Ganho pelo critério de módulo}")
        if pertence:
            doc_lines.append(eq(r"|K P(s_t)|=1\Rightarrow K=\frac{\prod_i|s_t-p_i|}{|C|\prod_j|s_t-z_j|}"))
            doc_lines.append(eq(rf"\boxed{{K={ptest['K']:.6f}}}"))
        else:
            doc_lines.append(r"Como a condição de ângulo não foi satisfeita, não se calcula $K$ para $K>0$.")

        # Gráfico
        try:
            graph_bytes = self._gerar_grafico_png_pdf()
            with tempfile.TemporaryDirectory() as td:
                td_path = Path(td)
                (td_path / "lgr.png").write_bytes(graph_bytes.getvalue())
                doc_lines.extend([
                    r"\passo{Diagrama final do LGR}",
                    r"\begin{center}",
                    r"\includegraphics[width=0.72\linewidth]{lgr.png}",
                    r"\end{center}",
                    r"\vfill",
                    r"\begin{center}\small\textit{Modo Prova: use as expressões acima como roteiro e reproduza no papel apenas o que o enunciado exigir.}\end{center}",
                    r"\end{document}",
                ])
                tex_path = td_path / "resolucao_LGR_modo_prova.tex"
                tex_path.write_text("\n".join(doc_lines), encoding="utf-8")
                result = subprocess.run(
                    [pdflatex, "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
                    cwd=td_path, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, timeout=30, check=False,
                )
                pdf_path = td_path / "resolucao_LGR_modo_prova.pdf"
                if result.returncode == 0 and pdf_path.exists():
                    return pdf_path.read_bytes()
                # Uma segunda compilação costuma resolver referências auxiliares.
                subprocess.run(
                    [pdflatex, "-interaction=nonstopmode", tex_path.name],
                    cwd=td_path, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, timeout=30, check=False,
                )
                if pdf_path.exists():
                    return pdf_path.read_bytes()
        except Exception:
            pass

        # Se o LaTeX estiver instalado mas a compilação falhar, não mascaramos
        # o problema devolvendo o PDF completo. O usuário ainda recebe um PDF
        # COMPACTO de Modo Prova, mantendo a finalidade do botão.
        return self._gerar_pdf_modo_prova_fallback(
            p1, p4, p7, p8, p9, p10, ptest, ponto_teste
        )

    def _gerar_pdf_modo_prova_fallback(
        self, p1, p4, p7, p8, p9, p10, ptest, ponto_teste,
    ) -> bytes:
        """Fallback compacto do Modo Prova quando pdflatex não está disponível.

        Esta versão NÃO chama ``gerar_pdf_resolucao``: mesmo sem LaTeX, o botão
        continua entregando um relatório enxuto, com as contas essenciais para
        reprodução manual.
        """
        font_regular, font_bold = self._pdf_font_setup()
        self._pdf_test_point = ponto_teste
        if self._last_rlist is None:
            self.calcular_lgr_exato()

        buf = BytesIO()
        doc = SimpleDocTemplate(
            buf, pagesize=A4,
            rightMargin=1.25 * cm, leftMargin=1.25 * cm,
            topMargin=1.15 * cm, bottomMargin=1.15 * cm,
            title="LGR - Modo Prova",
            author="Calculadora LGR",
        )
        styles = getSampleStyleSheet()
        title = ParagraphStyle(
            "MPTitle", parent=styles["Title"], fontName=font_bold,
            fontSize=16, leading=19, alignment=TA_CENTER, spaceAfter=7,
        )
        h1 = ParagraphStyle(
            "MPH1", parent=styles["Heading1"], fontName=font_bold,
            fontSize=10.5, leading=12.5, spaceBefore=5, spaceAfter=3,
        )
        body = ParagraphStyle(
            "MPBody", parent=styles["BodyText"], fontName=font_regular,
            fontSize=8.1, leading=10, spaceAfter=2.5,
        )
        eq = ParagraphStyle(
            "MPEq", parent=body, fontSize=8.3, leading=10.5,
            leftIndent=6, spaceAfter=2,
        )
        small = ParagraphStyle(
            "MPSmall", parent=body, fontSize=7.1, leading=8.4,
        )

        story = [
            Paragraph("Lugar Geométrico das Raízes — Modo Prova", title),
            Paragraph(
                "Versão compacta para reprodução manual. LaTeX indisponível: usando formato compacto de emergência.",
                small,
            ),
        ]

        def add_eq(text):
            story.append(Paragraph(escape(text), eq))

        # Dados + 1/2
        story.append(Paragraph("Dados", h1))
        add_eq(f"G(s) = {self._pdf_expr(self._poly_to_expr(self.numG))} / {self._pdf_expr(self._poly_to_expr(self.denG))}")
        add_eq(f"H(s) = {self._pdf_expr(self._poly_to_expr(self.numH))} / {self._pdf_expr(self._poly_to_expr(self.denH))}")
        add_eq(f"s_t = {self._fmt_complex(complex(ponto_teste), 4)}")

        story.append(Paragraph("1. Polinômio característico", h1))
        add_eq(f"1 + K·P(s) = 0")
        add_eq(f"Φ(s,K) = {self._pdf_expr(p1['char_expr'])} = 0")

        story.append(Paragraph("2. Pólos e zeros", h1))
        add_eq("Pólos: " + (", ".join(self._fmt_complex(p, 4) for p in self.polos) or "nenhum"))
        add_eq("Zeros: " + (", ".join(self._fmt_complex(z, 4) for z in self.zeros) or "nenhum"))
        add_eq(f"nP = {self.np}; nZ = {self.nz}")

        story.append(Paragraph("3–6. Plano-s, eixo real, ramos e simetria", h1))
        segs = p4.get("segmentos", [])
        seg_text = ", ".join(
            f"({'-∞' if math.isinf(a) else f'{a:.4f}'}, {'+∞' if math.isinf(b) else f'{b:.4f}'})"
            for a, b in segs
        ) or "nenhum"
        add_eq(f"Segmentos reais: {seg_text}")
        add_eq(f"n_ramos = nP = {self.np}; simétrico em relação ao eixo real.")

        story.append(Paragraph("7. Assíntotas", h1))
        if p7.get("numero", 0) > 0:
            angles = ", ".join(f"{a:.2f}°" for a in p7.get("angulos", []))
            add_eq(f"q = nP − nZ = {p7['numero']}; σA = {p7['centro']:.5f}; φA = {angles}")
        else:
            add_eq("Sem assíntotas para o infinito.")

        story.append(Paragraph("8. Saída/chegada", h1))
        add_eq("K(s) = −D(s)/N(s);   D'(s)N(s) − D(s)N'(s) = 0")
        valid = [c for c in p8.get("candidatos", []) if c.get("real") and c.get("valido")]
        if valid:
            add_eq("; ".join(
                f"s_b = {c['s'].real:.5f}, K_b = {c['K']:.6f}"
                for c in valid
            ))
        else:
            add_eq("Nenhum ponto válido para K > 0 nos segmentos reais.")

        story.append(Paragraph("9. Cruzamento do eixo imaginário — Routh", h1))
        table = p9.get("routh", {}).get("table", [])
        if table:
            max_cols = max(len(row) for _, row in table)
            data = [["Linha"] + [f"C{i+1}" for i in range(max_cols)]]
            for power, row in table:
                data.append([f"s^{int(power)}"] + [self._pdf_expr(v) for v in row] + ["0"]*(max_cols-len(row)))
            t = Table(data, repeatRows=1, hAlign="LEFT")
            t.setStyle(TableStyle([
                ("FONTNAME", (0,0), (-1,-1), font_regular),
                ("FONTNAME", (0,0), (-1,0), font_bold),
                ("FONTSIZE", (0,0), (-1,-1), 6.8),
                ("GRID", (0,0), (-1,-1), 0.3, colors.grey),
                ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
                ("LEFTPADDING", (0,0), (-1,-1), 3),
                ("RIGHTPADDING", (0,0), (-1,-1), 3),
            ]))
            story.append(t)
            story.append(Spacer(1, 3))
        crossings = p9.get("cruzamentos", [])
        if crossings:
            for c in crossings:
                add_eq(f"CRUZAMENTO: s = ±j{c['w']:.6f};   K = {c['K']:.6f}")
        else:
            add_eq("Nenhum cruzamento detectado para K > 0.")

        story.append(Paragraph("10. Ângulos de partida/chegada", h1))
        if p10.get("resultados"):
            for item in p10["resultados"]:
                tipo = "partida" if item["tipo"] == "partida" else "chegada"
                if item["tipo"] == "partida":
                    add_eq(
                        f"s = {self._fmt_complex(item['ponto'],4)}: θ_partida = 180° − ({item['soma_polos']:.3f}°) + ({item['soma_zeros']:.3f}°) = {item['angulo']:.3f}°"
                    )
                else:
                    add_eq(
                        f"s = {self._fmt_complex(item['ponto'],4)}: θ_chegada = 180° − ({item['soma_zeros']:.3f}°) + ({item['soma_polos']:.3f}°) = {item['angulo']:.3f}°"
                    )

        story.append(Paragraph("11. Condição de ângulo", h1))
        add_eq(f"∠P(s_t) = {ptest['fase_mod']:.4f}°  (alvo 180° mod 360°)")
        add_eq("RESULTADO: " + ("PERTENCE ao LGR para K > 0." if ptest["pertence"] else "NÃO PERTENCE ao LGR para K > 0."))

        story.append(Paragraph("12. Ganho K", h1))
        if ptest["pertence"]:
            add_eq(f"K = Π|s_t−p_i| / (|C|·Π|s_t−z_j|) = {ptest['K']:.6f}")
        else:
            add_eq("K não é calculado: condição de ângulo não satisfeita.")

        story.append(PageBreak())
        story.append(Paragraph("Diagrama final do LGR", h1))
        story.append(Image(self._gerar_grafico_png_pdf(), width=17.0*cm, height=12.0*cm))

        def footer(canvas, doc):
            canvas.saveState()
            canvas.setFont(font_regular, 7)
            canvas.drawString(1.25*cm, 0.65*cm, "Calculadora LGR — Modo Prova")
            canvas.drawRightString(A4[0]-1.25*cm, 0.65*cm, f"Página {doc.page}")
            canvas.restoreState()

        doc.build(story, onFirstPage=footer, onLaterPages=footer)
        return buf.getvalue()

    def _gerar_grafico_png_pdf(self):
        """Gera uma imagem estática do LGR para inserir no PDF."""
        if self._last_rlist is None:
            self.calcular_lgr_exato()

        fig, ax = plt.subplots(figsize=(7.1, 5.2), dpi=170)
        ax.axhline(0, linewidth=1.0)
        ax.axvline(0, linewidth=1.0)

        rlist = np.asarray(self._last_rlist)
        if rlist.ndim == 2:
            for i in range(rlist.shape[1]):
                ax.plot(np.real(rlist[:, i]), np.imag(rlist[:, i]), linewidth=1.7)

        # Segmentos reais
        for a, b in self.real_segments:
            xs = [a if not math.isinf(a) else b - self.span,
                  b if not math.isinf(b) else a + self.span]
            ax.plot(xs, [0, 0], linewidth=4)

        # Assíntotas
        raio = self.span * 1.5
        if self.assymptote_center is not None:
            for ang in self.assymptote_angles:
                rad = math.radians(ang)
                dx, dy = math.cos(rad), math.sin(rad)
                ax.plot(
                    [self.assymptote_center, self.assymptote_center + raio * dx],
                    [0, raio * dy],
                    linestyle="--", linewidth=1.1, alpha=0.65,
                )

        # Polos agrupados por multiplicidade.
        for group in self._group_points(self.polos):
            alpha = 0.42 if group["multiplicity"] > 1 else 1.0
            ax.scatter(
                [group["point"].real], [group["point"].imag],
                marker="x", s=95, linewidths=2.0, alpha=alpha,
            )
            if group["multiplicity"] > 1:
                ax.annotate(f"m={group['multiplicity']}",
                            (group["point"].real, group["point"].imag),
                            xytext=(5, 5), textcoords="offset points", fontsize=8)

        for group in self._group_points(self.zeros):
            ax.scatter(
                [group["point"].real], [group["point"].imag],
                marker="o", s=90, facecolors="none", linewidths=1.8,
            )

        # Pontos de saída/chegada
        for c in self.breakaway_candidates:
            if c.get("valido") and c.get("real"):
                ax.scatter([c["s"].real], [0], marker="s", s=40)

        # Cruzamentos imaginários
        for c in self.crossings:
            ax.scatter([0, 0], [c["w"], -c["w"]], marker="D", s=38)

        # Ponto de teste, quando disponível via último registro.
        if hasattr(self, "_pdf_test_point") and self._pdf_test_point is not None:
            pt = complex(self._pdf_test_point)
            ax.scatter([pt.real], [pt.imag], marker="*", s=135)

        span = max(self.span, 4.0)
        xr = (self.min_x - 0.15 * span, self.max_x + 0.15 * span)
        if xr[0] == xr[1]:
            xr = (-5, 5)
        yr_half = max(span * 0.62, 2.5)
        ax.set_xlim(*xr)
        ax.set_ylim(-yr_half, yr_half)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Eixo Real (Re)")
        ax.set_ylabel("Eixo Imaginário (Im)")
        ax.set_title("Lugar Geométrico das Raízes")
        ax.grid(True, alpha=0.18)
        fig.tight_layout()

        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf

    def gerar_pdf_resolucao(
        self, p1, p4, p7, p8, p9, p10, ptest, ponto_teste,
    ) -> bytes:
        """Gera um PDF completo da resolução calculada pelo aplicativo."""
        font_regular, font_bold = self._pdf_font_setup()
        self._pdf_test_point = ponto_teste
        if self._last_rlist is None:
            self.calcular_lgr_exato()

        buf = BytesIO()
        doc = SimpleDocTemplate(
            buf, pagesize=A4,
            rightMargin=1.5 * cm, leftMargin=1.5 * cm,
            topMargin=1.45 * cm, bottomMargin=1.45 * cm,
            title="Resolução de Lugar Geométrico das Raízes",
            author="Calculadora LGR",
        )
        styles = getSampleStyleSheet()
        title = ParagraphStyle(
            "LGRTitle", parent=styles["Title"], fontName=font_bold,
            fontSize=18, leading=22, alignment=TA_CENTER, spaceAfter=10,
        )
        h1 = ParagraphStyle(
            "LGRH1", parent=styles["Heading1"], fontName=font_bold,
            fontSize=13, leading=16, spaceBefore=8, spaceAfter=7,
        )
        h2 = ParagraphStyle(
            "LGRH2", parent=styles["Heading2"], fontName=font_bold,
            fontSize=10.5, leading=13, spaceBefore=5, spaceAfter=4,
        )
        body = ParagraphStyle(
            "LGRBody", parent=styles["BodyText"], fontName=font_regular,
            fontSize=8.8, leading=12, spaceAfter=4,
        )
        small = ParagraphStyle(
            "LGRSmall", parent=body, fontSize=7.3, leading=9.5,
        )
        eq = ParagraphStyle(
            "LGREq", parent=body, fontName=font_regular, fontSize=8.5,
            leading=11, leftIndent=10, spaceAfter=4,
        )

        story = []
        story.append(Paragraph("Resolução - Lugar Geométrico das Raízes (LGR)", title))
        story.append(Paragraph(
            "Documento gerado automaticamente a partir dos valores informados no aplicativo. "
            "A resolução apresenta as etapas intermediárias utilizadas para a reprodução manual.", body
        ))
        story.append(Spacer(1, 3))

        story.append(Paragraph("Dados do problema", h1))
        story.append(Paragraph(
            f"G(s) = {escape(self._pdf_expr(self._poly_to_expr(self.numG)))} / "
            f"({escape(self._pdf_expr(self._poly_to_expr(self.denG)))})", body
        ))
        story.append(Paragraph(
            f"H(s) = {escape(self._pdf_expr(self._poly_to_expr(self.numH)))} / "
            f"({escape(self._pdf_expr(self._poly_to_expr(self.denH)))})", body
        ))
        story.append(Paragraph(
            f"Ponto de teste: s_t = {escape(self._fmt_complex(complex(ponto_teste), 6))}", body
        ))

        # Passos 1-3
        story.append(Paragraph("Passo 1 - Polinômio característico", h1))
        story.append(Paragraph("1 + G(s)H(s) = 0", eq))
        story.append(Paragraph(
            f"1 + K·P(s) = 1 + K·({escape(self._pdf_expr(p1['gh']) )}) = 0", eq
        ))
        story.append(Paragraph(
            f"P(s) = {escape(self._pdf_expr(p1['num_expr']))} / {escape(self._pdf_expr(p1['den_expr']))}", body
        ))
        story.append(Paragraph(
            f"Forma fatorada: P(s) = C·Π(s-z_i) / Π(s-p_i), com C = {escape(self._pdf_expr(self.ganho_constante))}", body
        ))
        story.append(Paragraph(
            f"Φ(s,K) = {escape(self._pdf_expr(p1['char_expr']))} = 0", eq
        ))
        story.append(Paragraph(
            f"∂Φ/∂s = {escape(self._pdf_expr(p1['char_derivative']))}", eq
        ))

        story.append(Paragraph("Passo 2 - Pólos e zeros", h1))
        story.append(Paragraph(
            "Pólos: " + ", ".join(self._fmt_complex(p, 5) for p in self.polos) if self.polos else "Pólos: nenhum", body
        ))
        story.append(Paragraph(
            "Zeros: " + (", ".join(self._fmt_complex(z, 5) for z in self.zeros) if self.zeros else "nenhum"), body
        ))
        story.append(Paragraph(f"nP = {self.np}; nZ = {self.nz}.", body))

        story.append(Paragraph("Passo 3 - Esboço inicial no plano-s", h1))
        story.append(Paragraph("X = pólo; O = zero. O LGR inicia nos pólos e termina nos zeros (finitos ou infinitos).", body))

        # Passo 4
        story.append(Paragraph("Passo 4 - Segmentos do eixo real", h1))
        for item in p4.get("testes", []):
            story.append(Paragraph(
                f"Teste s = {item['ponto_teste']:.6f}: N_direita = {item['elementos_direita']} ({item['paridade']}). "
                f"Intervalo ({'-∞' if math.isinf(item['esquerda']) else f'{item['esquerda']:.6f}'}, "
                f"{'+∞' if math.isinf(item['direita']) else f'{item['direita']:.6f}'}) -> "
                f"{'pertence' if item['pertence'] else 'não pertence'} ao LGR.", small
            ))
        story.append(Paragraph(
            "Segmentos pertencentes: " + (", ".join(
                f"({'-∞' if math.isinf(a) else f'{a:.5f}'}, {'+∞' if math.isinf(b) else f'{b:.5f}'})"
                for a, b in p4.get("segmentos", [])
            ) if p4.get("segmentos") else "nenhum"), body
        ))

        story.append(Paragraph("Passo 5 - Número de lugares separados", h1))
        story.append(Paragraph(f"LS = nP = {self.np}. Logo, existem {self.np} ramos.", body))

        story.append(Paragraph("Passo 6 - Simetria", h1))
        story.append(Paragraph("Como os coeficientes do sistema são reais, o LGR é simétrico em relação ao eixo real.", body))

        # Passo 7
        story.append(Paragraph("Passo 7 - Assíntotas", h1))
        if p7.get("numero", 0) > 0:
            story.append(Paragraph(f"NA = nP - nZ = {p7['numero']}.", body))
            story.append(Paragraph(
                f"σA = (Σp_i - Σz_i)/(nP-nZ) = "
                f"({escape(self._pdf_expr(sp.expand(p7['soma_polos'])))} - {escape(self._pdf_expr(sp.expand(p7['soma_zeros'])))})/{p7['numero']} = {p7['centro']:.6f}",
                eq,
            ))
            for q, angle in enumerate(p7.get("angulos", [])):
                story.append(Paragraph(
                    f"φA,{q} = (2·{q}+1)·180°/{p7['numero']} = {angle:.6f}°", small
                ))
        else:
            story.append(Paragraph("nP - nZ ≤ 0: não há assíntotas para o infinito.", body))

        # Passo 8
        story.append(Paragraph("Passo 8 - Pontos de saída/chegada", h1))
        story.append(Paragraph(
            f"K(s) = -D(s)/N(s) = -({escape(self._pdf_expr(p8['K_num']))})/({escape(self._pdf_expr(p8['K_den']))})", eq
        ))
        story.append(Paragraph(f"dK/ds = {escape(self._pdf_expr(p8['dK_ds']))}", eq))
        story.append(Paragraph("dK/ds = 0 ⇔ D'(s)N(s) - D(s)N'(s) = 0", eq))
        story.append(Paragraph(
            f"D'(s) = {escape(self._pdf_expr(p8['Dp']))}; N'(s) = {escape(self._pdf_expr(p8['Np']))}.", body
        ))
        story.append(Paragraph(
            f"Equação numerador: {escape(self._pdf_expr(p8['equacao_numerador']))} = 0", eq
        ))
        for cand in p8.get("candidatos", []):
            if cand.get("real"):
                k = cand.get("K")
                ks = "indefinido" if k is None else f"{k:.6f}"
                status = "VÁLIDO" if cand.get("valido") else "ignorado"
                story.append(Paragraph(
                    f"Candidato s = {cand['s'].real:.6f}, K = {ks}: {status} para K>0 e segmento real.", small
                ))
            else:
                story.append(Paragraph(
                    f"Candidato complexo s = {escape(self._fmt_complex(cand['s'], 6))}: ignorado.", small
                ))

        # Passo 9
        story.append(Paragraph("Passo 9 - Cruzamento do eixo imaginário (Routh-Hurwitz)", h1))
        routh = p9.get("routh", {})
        table = routh.get("table", [])
        if table:
            max_cols = max(len(row) for _, row in table)
            data = [["Linha"] + [f"Coluna {i+1}" for i in range(max_cols)]]
            for power, row in table:
                data.append([f"s^{int(power)}"] + [self._pdf_expr(v) for v in row] + ["0"] * (max_cols - len(row)))
            t = Table(data, repeatRows=1, hAlign="LEFT")
            t.setStyle(TableStyle([
                ("FONTNAME", (0, 0), (-1, -1), "LGR-Regular"),
                ("FONTNAME", (0, 0), (-1, 0), "LGR-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 7),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.grey),
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ]))
            story.append(t)
            story.append(Spacer(1, 5))
        story.append(Paragraph(
            "Primeira coluna: " + ", ".join(
                f"s^{int(power)}: {escape(self._pdf_expr(row[0]))}" for power, row in table
            ), small
        ))
        story.append(Paragraph(
            "Valores positivos na primeira coluna determinam a faixa de estabilidade (para coeficiente líder positivo).", body
        ))
        if routh.get("K_candidates"):
            story.append(Paragraph(
                "Candidatos por anulação da primeira coluna: " + ", ".join(f"K = {k:.6f}" for k in routh["K_candidates"]), body
            ))
        if routh.get("used_epsilon"):
            story.append(Paragraph("Foi utilizado ε>0 devido a primeiro elemento nulo em uma linha da tabela de Routh.", small))

        s2 = p9.get("s2_method", {})
        if s2.get("disponivel"):
            story.append(Paragraph("Determinação do cruzamento pela equação da linha s²", h2))
            a, b = s2["linha_s2"]
            story.append(Paragraph(f"Linha s²: {escape(self._pdf_expr(a))}·s² + {escape(self._pdf_expr(b))} = 0", eq))
            story.append(Paragraph(f"Com s=jω: {escape(self._pdf_expr(s2['aux_jw']))}=0", eq))
            if s2.get("k_expression") is not None:
                story.append(Paragraph(f"Da equação de s²: K = {escape(self._pdf_expr(s2['k_expression']))}", eq))
            story.append(Paragraph(f"Re[Φ(jω,K)] = {escape(self._pdf_expr(s2['char_re']))}", eq))
            story.append(Paragraph(f"Im[Φ(jω,K)] = {escape(self._pdf_expr(s2['char_im']))}", eq))
            story.append(Paragraph(f"Im[Φ]/ω = 0 (ω≠0) ⇒ {escape(self._pdf_expr(s2['omega_equation']))} = 0", eq))
            for calc in s2.get("calculos", []):
                story.append(Paragraph(
                    f"ω = {calc['w']:.6f} ⇒ K = {calc['K']:.6f}; "
                    f"verificação Re[Φ] = {calc['real_residual']:.3e}, Im[Φ] = {calc['imag_residual']:.3e}.", small
                ))
        if p9.get("cruzamentos"):
            for c in p9["cruzamentos"]:
                story.append(Paragraph(
                    f"Cruzamento: s = ±j{c['w']:.6f}, K = {c['K']:.6f}.", body
                ))
        else:
            story.append(Paragraph("Não foi identificado cruzamento do eixo imaginário para K>0.", body))

        # Passo 10
        story.append(Paragraph("Passo 10 - Ângulos de partida e chegada", h1))
        if not p10.get("resultados"):
            story.append(Paragraph("Não há pólos/zeros complexos que exijam cálculo de ângulo.", body))
        for item in p10.get("resultados", []):
            tipo = "partida" if item["tipo"] == "partida" else "chegada"
            story.append(Paragraph(
                f"Ângulo de {tipo} em s = {escape(self._fmt_complex(item['ponto'], 5))}", h2
            ))
            vectors = item.get("vetores_polos", []) + item.get("vetores_zeros", [])
            if vectors:
                data = [["Vetor", "Origem", "Destino", "ΔRe", "ΔIm", "|v|", "Ângulo (°)"]]
                for i, v in enumerate(vectors, 1):
                    data.append([str(i), self._fmt_complex(v.origem, 4), self._fmt_complex(v.destino, 4),
                                 f"{v.dx:.4f}", f"{v.dy:.4f}", f"{v.magnitude:.4f}", f"{v.angulo_deg:.4f}"])
                vt = Table(data, repeatRows=1, hAlign="LEFT")
                vt.setStyle(TableStyle([
                    ("FONTNAME", (0, 0), (-1, -1), "LGR-Regular"),
                    ("FONTNAME", (0, 0), (-1, 0), "LGR-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 6.3),
                    ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                    ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ]))
                story.append(vt)
            if item["tipo"] == "partida":
                story.append(Paragraph(
                    f"Σθ_i = {item['soma_polos']:.4f}°; Σφ_j = {item['soma_zeros']:.4f}°; "
                    f"θ_partida = 180° - Σθ_i + Σφ_j = {item['angulo']:.4f}°.", body
                ))
            else:
                story.append(Paragraph(
                    f"Σφ_i = {item['soma_zeros']:.4f}°; Σθ_j = {item['soma_polos']:.4f}°; "
                    f"θ_chegada = 180° - Σφ_i + Σθ_j = {item['angulo']:.4f}°.", body
                ))

        # Passos 11 e 12
        story.append(Paragraph("Passo 11 - Teste da condição de ângulo", h1))
        story.append(Paragraph("Para cada vetor: θ = atan2(ΔIm, ΔRe) e |v| = √[(ΔRe)^2 + (ΔIm)^2].", body))
        vectors = ptest.get("vetores_polos", []) + ptest.get("vetores_zeros", [])
        if vectors:
            data = [["Vetor", "Origem", "Destino", "ΔRe", "ΔIm", "|v|", "Ângulo (°)"]]
            for i, v in enumerate(vectors, 1):
                data.append([str(i), self._fmt_complex(v.origem, 4), self._fmt_complex(v.destino, 4),
                             f"{v.dx:.4f}", f"{v.dy:.4f}", f"{v.magnitude:.4f}", f"{v.angulo_deg:.4f}"])
            vt = Table(data, repeatRows=1, hAlign="LEFT")
            vt.setStyle(TableStyle([
                ("FONTNAME", (0, 0), (-1, -1), "LGR-Regular"),
                ("FONTNAME", (0, 0), (-1, 0), "LGR-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 6.3),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
            ]))
            story.append(vt)
        story.append(Paragraph(
            f"arg(C) = {ptest['fase_constante_deg']:.4f}°; Σθ_i = {ptest['soma_angulos_polos']:.4f}°; "
            f"Σφ_j = {ptest['soma_angulos_zeros']:.4f}°.", body
        ))
        story.append(Paragraph(
            f"∠P(s_t) = arg(C) + Σφ_j - Σθ_i = {ptest['fase_bruta']:.4f}° ≡ {ptest['fase_mod']:.4f}° (mod 360°).", eq
        ))
        story.append(Paragraph(
            "Resultado: ponto pertence ao LGR para K>0." if ptest["pertence"] else "Resultado: ponto não pertence ao LGR para K>0.", body
        ))

        story.append(Paragraph("Passo 12 - Determinação do ganho K pela condição de módulo", h1))
        if ptest["pertence"]:
            story.append(Paragraph(
                "|K·P(s_t)| = 1 ⇒ K = Π|s_t-p_i| / (|C|·Π|s_t-z_j|).", eq
            ))
            story.append(Paragraph(
                f"K = {ptest['produto_mod_polos']:.6f} / "
                f"(|{escape(self._pdf_expr(ptest['ganho_constante']))}|·{ptest['produto_mod_zeros']:.6f}) "
                f"= {ptest['K']:.6f}.", eq
            ))
            story.append(Paragraph(
                f"Logo, para s_t = {escape(self._fmt_complex(ponto_teste, 6))}, K = {ptest['K']:.6f}.", body
            ))
        else:
            story.append(Paragraph("Como a condição de ângulo não foi satisfeita, não se calcula K para K>0.", body))

        # Gráfico final
        story.append(Paragraph("Diagrama final do LGR", h1))
        story.append(Image(self._gerar_grafico_png_pdf(), width=17.2*cm, height=12.1*cm))

        story.append(Spacer(1, 8))
        story.append(Paragraph(
            "Fim da resolução. Os valores e cálculos deste documento correspondem aos dados informados no aplicativo.", small
        ))

        def footer(canvas, doc):
            canvas.saveState()
            canvas.setFont("LGR-Regular", 7)
            canvas.drawString(1.5*cm, 0.75*cm, "Calculadora LGR - resolução automática")
            canvas.drawRightString(A4[0]-1.5*cm, 0.75*cm, f"Página {doc.page}")
            canvas.restoreState()

        doc.build(story, onFirstPage=footer, onLaterPages=footer)
        return buf.getvalue()

    def _poly_to_expr(self, coeffs):
        """Converte lista de coeficientes em polinômio simbólico."""
        return sp.Poly.from_list(list(np.asarray(coeffs, dtype=float)), gens=self.s).as_expr()
