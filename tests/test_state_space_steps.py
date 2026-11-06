from types import SimpleNamespace
import sympy as sp
import pytest
from core.state_space_steps import analysis_steps, placement_steps, canonical_steps, build_resolution, resolution_markdown
S=sp.Symbol('s')

@pytest.fixture()
def system():
    A=sp.Matrix([[-4,3,0],[-1,-4,-1],[0,0,-2]])
    B=sp.Matrix([0,0,3]); C=sp.Matrix([[1,0,0]])
    U=sp.Matrix.hstack(B,A*B,A*A*B)
    V=sp.Matrix.vstack(C,C*A,C*A*A)
    a=SimpleNamespace(A=A,B=B,C=C,U=U,V=V,rank_U=U.rank(),rank_V=V.rank(),characteristic=A.charpoly(S).as_expr(),eigenvalues=(-4-3**.5*1j,-4+3**.5*1j,-2+0j),stable=True,controllable=True,observable=True)
    wanted=sp.expand((S+10)*(S*S+4*S+5)); qA=A*A*A+14*A*A+45*A+50*sp.eye(3)
    e=sp.Matrix([[0,0,1]]);K=-e*U.inv()*qA
    Ac=A+B*K
    p=SimpleNamespace(mode='K',desired=wanted,qA=qA,gain=K,result_A=Ac,result_characteristic=Ac.charpoly(S).as_expr(),result_eigenvalues=(-10+0j,-2+1j,-2-1j),inverse_matrix=U.inv())
    return a,p

def test_analysis_demonstrates_determinant_expansion_and_matrix_products(system):
    a,p=system
    steps=analysis_steps(a)
    tex=[x.value for x in steps if x.kind=='math']
    assert any('\\det(sI-A)' in s for s in tex)
    assert sum('\\Delta(s)' in s for s in tex)>=3
    assert any('AB=' in s for s in tex)
    assert any('CA=' in s for s in tex)
    assert any('\\det(U)' in s for s in tex)
    assert any('\\det(V)' in s for s in tex)

def test_feedback_reconstructible_numeric_substitution(system):
    a,p=system
    steps=placement_steps(a,p)
    equations='\n'.join(x.value for x in steps if x.kind=='math')
    assert 'q(A)=A^{3}+14A^{2}+45A+50I' in equations
    assert 'U^{-1}=' in equations
    assert 'e_n^TU^{-1}=' in equations
    assert 'A_{mf}=A+BK=' in equations
    assert r'\frac{8}{3}' in equations
    assert r'0&\cdots&0' not in equations # vector de seleção explícito

def test_observer_stages(system):
    a,p=system
    targets=sp.expand((S+10)*(S*S+2*S+2));q=A_mat(a)**3 + 12*A_mat(a)**2+22*A_mat(a)+20*sp.eye(3)
    L=q*a.V.inv()*sp.Matrix([0,0,1]);Aobs=a.A-L*a.C
    ob=SimpleNamespace(mode='L',desired=targets,qA=q,gain=L,result_A=Aobs,result_characteristic=Aobs.charpoly(S).as_expr(),result_eigenvalues=(-10+0j,-1+1j,-1-1j),inverse_matrix=a.V.inv())
    steps=placement_steps(a,ob)
    eq='\n'.join(x.value for x in steps if x.kind=='math')
    assert r'L=q(A)V^{-1}e_n' in eq
    assert 'V^{-1}e_n=' in eq
    assert 'A_o=A-LC=' in eq

def A_mat(a): return a.A

def test_input_shown_once_and_export_matches(system):
    a,p=system
    steps=build_resolution(a,p,operation='K')
    eq='\n'.join(x.value for x in steps if x.kind=='math')
    assert eq.count('A='+sp.latex(a.A))==1
    md=resolution_markdown(steps)
    assert '\\Delta_d(s)' in md and 'e_n^TU^{-1}' in md
    assert md.count('$$') >= 4

def test_canonical_coefficient_synthesis():
    A=sp.Matrix([[0,1,0],[0,0,1],[-5,-7,-3]])
    B=sp.Matrix([0,0,1]); C=sp.Matrix([[-2,4,3]])
    g=(3*S*S+4*S-2)/(S**3+3*S*S+7*S+5)
    form=SimpleNamespace(A=A,B=B,C=C,D=0,original_transfer=g,reconstructed_transfer=g)
    steps=canonical_steps(form)
    eq='\n'.join(x.value for x in steps if x.kind=='math')
    assert r'\alpha_{1}=3' in eq
    assert r'\beta_{3}=-2' in eq
    assert 'A_c=' in eq and 'B_c=' in eq and 'C_c=' in eq
    assert 'G(s)' in eq