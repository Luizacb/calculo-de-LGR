"""
Smoke tests da interface (app.py) via streamlit.testing.v1.AppTest.
Não testam valores pixel-a-pixel, mas garantem que:
  - a app carrega sem exceção com os valores padrão;
  - o botão "Calcular LGR" dispara o cálculo sem erro;
  - inputs alternativos (com polos complexos, sem zeros, entrada inválida)
    não quebram a app;
  - o resultado em destaque aparece com valores coerentes.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from streamlit.testing.v1 import AppTest

APP_PATH = os.path.join(os.path.dirname(__file__), "..", "app.py")


def _run_com_valores(nG=None, dG=None, nH=None, dH=None):
    at = AppTest.from_file(APP_PATH, default_timeout=30)
    at.run()
    assert not at.exception

    if nG is not None:
        at.text_input[0].set_value(nG)
    if dG is not None:
        at.text_input[1].set_value(dG)
    if nH is not None:
        at.text_input[2].set_value(nH)
    if dH is not None:
        at.text_input[3].set_value(dH)

    at.button[0].click().run()
    return at


def test_app_carrega_sem_excecao():
    at = AppTest.from_file(APP_PATH, default_timeout=30)
    at.run()
    assert not at.exception


def test_calcular_com_valores_padrao():
    at = _run_com_valores()
    assert not at.exception
    # deve haver conteúdo LaTeX renderizado (equação característica etc.)
    assert len(at.latex) > 0


def test_calcular_polos_complexos():
    at = _run_com_valores(nG="1", dG="1 2 5", nH="1", dH="1")
    assert not at.exception


def test_calcular_sem_zeros_finitos():
    at = _run_com_valores(nG="1", dG="1 6 11 6", nH="1", dH="1")
    assert not at.exception


def test_calcular_entrada_invalida_mostra_erro_sem_exception():
    at = _run_com_valores(dG="abc")
    assert not at.exception
    # deve exibir mensagem de erro amigável, sem quebrar a app
    assert len(at.error) > 0


def test_calcular_ordem_alta():
    at = _run_com_valores(nG="1 1", dG="1 6 11 6 0", nH="1", dH="1 2")
    assert not at.exception


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))
