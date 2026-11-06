"""Cobertura de presets completos e relatórios PDF/HTML."""
from io import BytesIO
from zipfile import ZipFile

import pytest

from core.state_space_questions import QUESTION_TITLES, solve_full_question, individual_document
from core.state_space_steps import ResolutionStep
from reports.state_space_export import (build_state_space_html, build_state_space_pdf,
                                        build_all_state_space_zip)


@pytest.mark.parametrize('question,sections', [(1,2),(2,3),(3,3),(4,3)])
def test_question_is_complete_and_has_ordered_sections(question,sections):
    doc = solve_full_question(question)
    assert doc.title == QUESTION_TITLES[question]
    assert len(doc.sections) == sections
    text = '\n'.join(s.value for part in doc.sections for s in part.steps)
    assert 'Delta' in text or r'\Delta' in text
    if question == 1:
        assert 'Sistema 1' in doc.sections[0].title
        assert 'Sistema 2' in doc.sections[1].title
        assert sum(step.kind == "heading" and step.value == "Conferência dos dados de entrada"
                   for sec in doc.sections for step in sec.steps) == 2
    if question in (2,3):
        assert sum('Conferência dos dados' in s.value for section in doc.sections for s in section.steps)==1
        assert 'K=' in text and 'L=' in text
    if question == 4:
        assert 'Forma canônica' in doc.sections[0].title
        assert 'Observador' in doc.sections[1].title
        assert 'Realimentação' in doc.sections[2].title


def test_pdf_html_exports_match_sections_and_are_renderable():
    doc = solve_full_question(1)
    pdf = build_state_space_pdf(doc)
    html = build_state_space_html(doc).decode('utf-8')
    assert pdf.startswith(b'%PDF-') and len(pdf)>1000
    assert '<!doctype html>' in html and 'MathJax' in html
    assert 'Sistema 1' in html and 'Sistema 2' in html
    assert 'Controlabilidade' in html and 'Observabilidade' in html


def test_individual_mode_exports_only_one_part():
    doc = individual_document('Exercício isolado', [ResolutionStep('math',r'x=1')])
    assert len(doc.sections)==1
    assert b'%PDF-' in build_state_space_pdf(doc)[:8]


def test_zip_has_exactly_four_pdfs_and_four_htmls():
    with ZipFile(BytesIO(build_all_state_space_zip(solve_full_question))) as z:
        names=z.namelist()
        assert len(names)==8
        assert sorted(names)==sorted(f'U3L1_Questao_{k}.{ext}' for k in range(1,5) for ext in ('pdf','html'))
        for k in range(1,5):
            assert z.read(f'U3L1_Questao_{k}.pdf').startswith(b'%PDF-')
            assert b'MathJax' in z.read(f'U3L1_Questao_{k}.html')