"""PDF e HTML completos da Unidade 3, Lista 1.

Sem depender de executáveis LaTeX ou de serviços externos para gerar o PDF.
A sequência matemática é compartilhada com a interface, evitando divergências.
"""
from html import escape
from io import BytesIO
import re
from zipfile import ZIP_DEFLATED, ZipFile

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import Image, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


_STYLE = """
:root {color-scheme:light} * {box-sizing:border-box}
html,body {background:#fff!important;color:#172333!important;font-family:Arial,Helvetica,sans-serif}
body {max-width:1040px;margin:auto;padding:38px 46px;line-height:1.55}
header {border-bottom:3px solid #244574;margin-bottom:28px;padding-bottom:16px}
header small {display:block;letter-spacing:.08em;color:#53667f;text-transform:uppercase}
h1 {font-size:27px;margin:6px 0 9px} h2 {font-size:21px;margin:32px 0 13px;color:#183c66}
h3 {font-size:16px;margin:23px 0 10px} .step {margin:13px 0;}
.math {margin:12px 0;overflow-x:auto;padding:4px 8px;}
.note {color:#506478;font-size:14px}.success {border-left:3px solid #367b58;padding:7px 12px}
.warning {border-left:3px solid #bf682d;padding:7px 12px}.detail {color:#52657b}
section.question {padding:4px 0 14px;border-bottom:1px solid #dce4ed}
footer {margin-top:34px;color:#65768d;font-size:12px}
@media print {body {max-width:none;padding:0;color:#000} section.question {break-before:auto;border:none}
h2,h3 {break-after:avoid} .math {break-inside:avoid}}
"""


def _iter_steps(document):
    for section in document.sections:
        yield section.title, section.steps


def build_state_space_html(document):
    """HTML autônomo quanto ao texto; MathJax via CDN para fórmulas visualmente completas."""
    label = getattr(document, "list_label", "U3L1")
    list_name = "Lista 2" if label == "U3L2" else "Lista 1"
    parts = [
        '<!doctype html><html lang="pt-BR"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width,initial-scale=1">',
        '<title>' + escape(document.title) + '</title>',
        '<style>' + _STYLE + '</style>',
        '<script>window.MathJax={tex:{displayMath:[["\\\\[","\\\\]"]]},svg:{fontCache:"global"}};</script>',
        '<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>',
        '</head><body><header><small>Projeto de Sistemas de Controle · UFRN · 2026.2</small>',
        '<h1>' + escape(document.title) + '</h1>',
        '<p class="note">Unidade 3, ' + list_name + ' · Método da apostila, capítulo 7 · ',
        'Realimentação u(t)=Kx(t), A+BK; observador A−LC.</p></header>',
    ]
    for section_title, steps in _iter_steps(document):
        parts.append('<section class="question"><h2>' + escape(section_title) + '</h2>')
        for step in steps:
            if step.kind == 'heading':
                parts.append('<h3>' + escape(step.value) + '</h3>')
            elif step.kind in ('math', 'detail_math'):
                classes = 'math detail' if step.kind == 'detail_math' else 'math'
                # Escaping HTML while leaving original TeX inside the MathJax delimiters.
                parts.append(f'<div class="{classes}">\\[{escape(step.value)}\\]</div>')
            else:
                css = step.kind if step.kind in ('note','success','warning') else 'note'
                parts.append(f'<p class="{css}">{escape(step.value)}</p>')
        parts.append('</section>')
    parts.append('<footer>Relatório gerado pela Calculadora de Controle · ' + label + ' · ',)
    parts.append('Equações verificadas simbolicamente com SymPy.</footer></body></html>')
    return ''.join(parts).encode('utf-8')


def _font_names():
    """Evita dependência obrigatória de fontes externas no Streamlit Cloud."""
    import os
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    regular_candidates = [
        '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
        '/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf',
    ]
    mono_candidates = [
        '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf',
        '/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf',
    ]
    family = []
    for font, candidates, fallback in [('U3Sans', regular_candidates, 'Helvetica'),
                                       ('U3Mono', mono_candidates, 'Courier')]:
        path = next((x for x in candidates if os.path.isfile(x)), None)
        if path:
            if font not in pdfmetrics.getRegisteredFontNames():
                pdfmetrics.registerFont(TTFont(font, path))
            family.append(font)
        else:
            family.append(fallback)
    return family


def _group_to_plain(src):
    """Frações simples e matrizes TeX em texto legível, sem perder coeficientes."""
    src = src.replace(r'\left', '').replace(r'\right', '')
    # Matriz convertida para linhas com separadores explícitos.
    pattern = re.compile(r'\\begin\{(?:bmatrix|matrix|pmatrix)\}(.*?)\\end\{(?:bmatrix|matrix|pmatrix)\}', re.S)
    def replace_matrix(match):
        rows = [row.strip() for row in match.group(1).split(r'\\') if row.strip()]
        return '[ ' + ' ; '.join(' , '.join(x.strip() for x in row.split('&')) for row in rows) + ' ]'
    src = pattern.sub(replace_matrix, src)
    # Frações aninhadas de baixa profundidade: converter das mais internas.
    for _ in range(8):
        changed = re.sub(r'\\frac\{([^{}]*)\}\{([^{}]*)\}',r'(\1)/(\2)',src)
        if changed == src:
            break
        src = changed
    replacements = {
        r'\Delta': 'Delta', r'\lambda': 'lambda', r'\alpha': 'alpha',
        r'\beta': 'beta', r'\omega': 'omega', r'\theta': 'theta',
        r'\sigma': 'sigma', r'\phi':'phi', r'\det': 'det',
        r'\operatorname{posto}':'posto',r'\operatorname':'',
        r'\cdot':'*',r'\quad':'  ',r'\qquad':'  ',r'\;':' ',
        r'\Longrightarrow': ' => ',r'\Rightarrow': ' => ',r'\Longleftrightarrow':' <=> ',
        r'\neq':' != ',r'\approx':' ~= ',r'\pm':' +/- ',r'\times':' x ',
        r'\cdots':'...',r'\vdots':'...',r'\mathrm':'',r'\text':'',
        r'\prod':'produto',r'\in':' pertence a ',r'\sqrt':'raiz',
        r'\{':'{',r'\}':'}',r'\,':' ',
    }
    for before, after in replacements.items():
        src = src.replace(before, after)
    src = re.sub(r'\\([A-Za-z]+)', r'\1', src)  # comando não reconhecido
    src = re.sub(r'\^\{([^{}]*)\}',r'^\1',src)
    src = re.sub(r'_\{([^{}]*)\}',r'_\1',src)
    src = src.replace('{','').replace('}','')
    src = src.replace('−','-')
    return ' '.join(src.split())


def _math_figure(tex, max_width):
    """Tiposet opcional para fórmulas sem matrizes; fallback é texto acessível."""
    if r'\begin{' in tex or len(tex)>220 or r'\operatorname' in tex:
        return None
    try:
        from matplotlib.mathtext import math_to_image
        buffer = BytesIO()
        math_to_image('$'+tex+'$', buffer, dpi=170, format='png', color='#202c42')
        buffer.seek(0)
        from PIL import Image as PILImage
        with PILImage.open(buffer) as img:
            w,h = img.size
        buffer.seek(0)
        natural_width = w*72/170
        scale = min(1.0, max_width/max(natural_width,1))
        return Image(buffer, width=natural_width*scale, height=h*72/170*scale)
    except Exception:
        return None



def _matrix_flowables(tex, math_style, mono_font):
    """Mostra cada matriz em linhas/colunas em vez de serializar todas as entradas."""
    matrix_regex = re.compile(r'\\begin\{(?:bmatrix|matrix|pmatrix)\}(.*?)\\end\{(?:bmatrix|matrix|pmatrix)\}', re.S)
    matches=list(matrix_regex.finditer(tex))
    if not matches:
        return []
    output=[]
    offset=0
    for match in matches:
        label=_group_to_plain(tex[offset:match.start()]).strip(' ,')
        label=re.sub(r'[\[(\s]+$', '', label)
        label=re.sub(r'^[\])\s,]+', '', label)
        if label:
            output.append(Paragraph(escape(label),math_style))
        rows = [[_group_to_plain(value.strip()) for value in row.split('&')]
                for row in match.group(1).split(r'\\') if row.strip()]
        if rows:
            table=Table(rows,hAlign='LEFT')
            table.setStyle(TableStyle([
                ('FONTNAME',(0,0),(-1,-1),mono_font),
                ('FONTSIZE',(0,0),(-1,-1),9),
                ('LEADING',(0,0),(-1,-1),13),
                ('ALIGN',(0,0),(-1,-1),'CENTER'),
                ('BACKGROUND',(0,0),(-1,-1),colors.HexColor('#f5f7fa')),
                ('LEFTPADDING',(0,0),(-1,-1),12),
                ('RIGHTPADDING',(0,0),(-1,-1),12),
                ('TOPPADDING',(0,0),(-1,-1),4),
                ('BOTTOMPADDING',(0,0),(-1,-1),4),
                ('LINEBEFORE',(0,0),(0,-1),0.85,colors.HexColor('#64748b')),
                ('LINEAFTER',(-1,0),(-1,-1),0.85,colors.HexColor('#64748b')),
            ]))
            output.append(table)
            output.append(Spacer(1,4))
        offset=match.end()
    trailing=_group_to_plain(tex[offset:]).strip(' ,')
    if trailing and trailing not in (')',']',')='):
        output.append(Paragraph(escape(trailing),math_style))
    return output

def _page_number(canvas, doc, label='U3L1'):
    canvas.saveState()
    canvas.setFont('Helvetica', 8)
    canvas.setFillColor(colors.HexColor('#5d7089'))
    canvas.drawString(1.7*cm, 1.0*cm, 'UFRN · Projeto de Sistemas de Controle · '+label)
    canvas.drawRightString(A4[0]-1.7*cm, 1.0*cm, f'{doc.page}')
    canvas.restoreState()


def build_state_space_pdf(document):
    """PDF multipágina com todas as contas e dados, sem dependências de LaTeX."""
    label = getattr(document, "list_label", "U3L1")
    sans, mono = _font_names()
    base = getSampleStyleSheet()
    styles = {
        'title': ParagraphStyle('U3Title', parent=base['Title'], fontName=sans,
                                fontSize=16, leading=21, spaceAfter=11,textColor=colors.HexColor('#15385d')),
        'section': ParagraphStyle('U3Section', parent=base['Heading2'], fontName=sans,
                                fontSize=12, leading=16, spaceBefore=16,spaceAfter=6,textColor=colors.HexColor('#244574')),
        'heading': ParagraphStyle('U3Heading', parent=base['Heading3'],fontName=sans,
                                fontSize=10.5,leading=13,spaceBefore=11,spaceAfter=5),
        'body': ParagraphStyle('U3Body',parent=base['BodyText'],fontName=sans,
                                fontSize=9.5,leading=13,spaceAfter=6),
        'math': ParagraphStyle('U3Math',parent=base['BodyText'],fontName=mono,
                                fontSize=8.0,leading=11,spaceAfter=8,wordWrap='CJK',
                                borderPadding=(3,4,3,4),backColor=colors.HexColor('#f5f7fa')),
        'foot': ParagraphStyle('U3Foot',parent=base['BodyText'],fontName=sans,
                                fontSize=8,leading=10,textColor=colors.grey),
    }
    buffer=BytesIO()
    doc = SimpleDocTemplate(buffer,pagesize=A4,
                             leftMargin=1.7*cm,rightMargin=1.7*cm,
                             topMargin=1.7*cm,bottomMargin=1.7*cm,
                             title=document.title,author='Calculadora de Controle '+label)
    story=[Paragraph(escape(document.title),styles['title']),
           Paragraph('UFRN · Projeto de Sistemas de Controle · 2026.2',styles['foot']),
           Paragraph('Apostila, capítulo 7. Convenção: u(t)=Kx(t) ⇒ A+BK; observador A-LC.',styles['body'])]
    available = A4[0]-3.4*cm-10
    for section_title, steps in _iter_steps(document):
        story.append(Paragraph(escape(section_title),styles['section']))
        for step in steps:
            if step.kind=='heading':
                story.append(Paragraph(escape(step.value),styles['heading']))
            elif step.kind in ('math','detail_math'):
                # A forma tipográfica é usada sempre que o motor local a suporta.
                if r'\begin{' in step.value:
                    story.extend(_matrix_flowables(step.value,styles['math'],mono))
                else:
                    graphic = _math_figure(step.value, available)
                    if graphic is not None:
                        story.append(graphic)
                        story.append(Spacer(1,4))
                    else:
                        story.append(Paragraph(escape(_group_to_plain(step.value)),styles['math']))
            else:
                story.append(Paragraph(escape(step.value),styles['body']))
    footer=lambda canvas, page: _page_number(canvas, page, label)
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    return buffer.getvalue()


def build_all_state_space_zip(solve_fn):
    """Pacote pronto para entrega: exatamente quatro PDFs e quatro HTMLs."""
    buffer=BytesIO()
    with ZipFile(buffer,'w',compression=ZIP_DEFLATED) as zipfile:
        for number in range(1,5):
            document=solve_fn(number)
            name=f'{getattr(document, "list_label", "U3L1")}_Questao_{number}'
            zipfile.writestr(name+'.pdf',build_state_space_pdf(document))
            zipfile.writestr(name+'.html',build_state_space_html(document))
    return buffer.getvalue()