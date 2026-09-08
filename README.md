# 🌸 LGR — Calculadora de Lugar Geométrico das Raízes

Aplicativo web de apoio para exercícios/provas de Sistemas de Controle
(cálculo do LGR — Lugar Geométrico das Raízes).

## Estrutura do projeto

```
lgr_rosa/
├── app.py                      # Interface Streamlit — só exibe, não calcula
├── core/
│   ├── __init__.py
│   └── calculos.py             # TODA a matemática do LGR (sem depender de Streamlit)
├── reference/
│   └── original_calculos.py    # Funções puras extraídas do repositório original,
│                                # usadas apenas como "gabarito" nos testes de regressão
├── tests/
│   ├── test_regression.py      # 124 testes: core.calculos  vs  reference.original_calculos
│   └── test_app_smoke.py       # 6 testes: a interface roda sem exceção (via AppTest)
├── requirements.txt
├── requirements-dev.txt        # requirements.txt + pytest
└── README.md
```

`core/calculos.py` e `reference/original_calculos.py` são deliberadamente
separados: o primeiro é a implementação nova (documentada, com type hints,
reorganizada), o segundo é o código original **inalterado**, mantido só
para o teste automatizado comparar os dois a qualquer momento que o
código evoluir.

## Rodando localmente

```bash
git clone <seu-repositório>
cd lgr_rosa
python3 -m venv .venv && source .venv/bin/activate   # opcional, recomendado
pip install -r requirements.txt
streamlit run app.py
```

O app abre em `http://localhost:8501`.

### Rodando os testes

```bash
pip install -r requirements-dev.txt
pytest tests/ -v
```

## Deploy no Streamlit Community Cloud

1. Suba este diretório para um repositório no GitHub (público ou privado).
2. Acesse [share.streamlit.io](https://share.streamlit.io), clique em
   **New app**.
3. Selecione o repositório, a branch e defina `app.py` como arquivo
   principal.
4. Clique em **Deploy** — o Streamlit Cloud instala automaticamente as
   dependências listadas em `requirements.txt`.

## Funcionalidades reproduzidas

Tudo que o app calcula e exibe:

- Montagem da equação característica `D(s) + K·N(s) = 0` a partir de
  G(s) = K·N_G/D_G e H(s) = N_H/D_H
- Polos e zeros de malha aberta (`numpy.roots`)
- Forma fatorada de P(s)
- Segmentos do eixo real pertencentes ao LGR
- Número de lugares separados e simetria
- Assíntotas (número, centroide σₐ, ângulos)
- Pontos de breakaway/break-in (via `dK/ds=0`), com classificação de
  pertinência ao LGR
- Tabela de Routh-Hurwitz simbólica (K como parâmetro), condições de
  estabilidade e valores críticos de K
- Cruzamento com o eixo imaginário (`s=jω`), com K e ω correspondentes
- Ângulos de partida (polos complexos) e de chegada (zeros complexos)
- Critério do ângulo para um ponto de teste s₀ (pertence/não pertence ao
  LGR, com a mesma tolerância de 5° do original)
- Critério do módulo — cálculo de K no ponto de teste
- Todos os 7 gráficos individuais + o gráfico consolidado do LGR completo
- Limites de gráfico manuais (opcional)
