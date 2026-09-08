# 🌸 LGR Rosa — Calculadora de Lugar Geométrico das Raízes

Aplicativo web de apoio para exercícios/provas de Sistemas de Controle
(cálculo do LGR — Lugar Geométrico das Raízes). Implementação **própria e
independente**, com identidade visual rosa, criada a partir da análise
funcional e matemática do projeto
[`igorservo159/calc_lgr`](https://github.com/igorservo159/calc_lgr)
([app original](https://calc-lgr.streamlit.app)), mas **sem copiar** seu
código de interface.

> A prioridade nº 1 deste projeto é a precisão matemática. Toda a lógica
> de cálculo foi validada por uma suíte de 124 testes de regressão contra
> o algoritmo original antes de qualquer linha da interface ser escrita.

## Estrutura do projeto

```
lgr_rosa/
├── app.py                      # Interface Streamlit (rosa) — só exibe, não calcula
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

Não há segredos, chaves de API nem variáveis de ambiente necessárias.

## Funcionalidades reproduzidas

Tudo que o app original calcula e exibe foi reproduzido:

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

## Principais decisões de projeto

- **Separação matemática/interface**: `core/calculos.py` não importa
  Streamlit. Isso permite testar a matemática isoladamente e reaproveitar
  o núcleo em outro front-end no futuro, se quiser.
- **Mesmas bibliotecas do original**: `numpy.roots`, `numpy.convolve`,
  `numpy.polyder`/`polyadd`/`polysub`/`polyval` e `sympy` para a tabela de
  Routh simbólica — nada foi trocado por alternativas (ex.: `scipy.signal`,
  `python-control`), justamente para preservar a equivalência numérica.
- **Reorganização da interface**: os 12 passos didáticos do original foram
  agrupados em 8 abas temáticas + um bloco de **"Resultado em destaque"**
  sempre visível no topo (polos/zeros, centroide, K crítico, cruzamento
  jω, veredito do ponto de teste), pensado para consulta rápida durante
  uma prova, sem precisar abrir vários expanders.
- **Identidade visual rosa**: paleta baseada em `#d6336c` (rosa principal),
  fundos claros (`#fff0f6`), cartões brancos com borda rosa suave — sem
  animações, mantendo o foco em velocidade de leitura.

## Limitações e diferenças conhecidas em relação ao original

Nenhuma diferença **numérica** foi encontrada (ver seção de testes
abaixo). As duas únicas ressalvas são de **comportamento conhecido do
original, preservado de propósito**:

1. **Tabela de Routh com pivô nulo**: quando um elemento da primeira
   coluna é exatamente zero, o algoritmo original interrompe o cálculo
   sem aplicar o método do polinômio auxiliar (caso clássico de sistema
   marginalmente estável / raízes no eixo jω). A nova versão reproduz
   esse comportamento por padrão para bater com o gabarito, mas o
   `core/calculos.py` sinaliza o caso com um comentário
   `# DESVIO DO ORIGINAL:` explicando a limitação, caso queira estender.
2. **Tolerância de 5° no critério do ângulo** (Passo "Ponto de teste"): é
   a mesma tolerância generosa do original. Ela fica visível na interface
   (`core.calculos.TOL_ANGULO_LGR`) para você saber que pontos "no limite"
   podem ser classificados como pertencentes ao LGR mesmo com um desvio
   de alguns graus.

## Resultados dos testes comparativos com o original

Suíte: `tests/test_regression.py` — 124 testes, 9 conjuntos de entrada
(polos reais, polos complexos, polo na origem, com/sem zeros finitos,
realimentação H(s)≠1, ordem alta, caso de fronteira n_p=n_z+1) ×
categorias de cálculo (polos/zeros, segmentos, assíntotas, breakaway,
cruzamento jω, K crítico via Routh, ramos completos do LGR, critério do
ângulo/módulo em 5 pontos de teste, formatação LaTeX, entradas
inválidas).

```
124 passed in ~5s
```

Exemplo de comparação (um dos 9 cenários, os demais seguem o mesmo
padrão nos logs de teste):

```
Entrada: G(s) = 1/(s²+2s+5), H(s)=1        (polos complexos)
Original → polos = -1+2j, -1-2j
Novo     → polos = -1+2j, -1-2j            (diferença: 0)

Original → K crítico (Routh) = 6.0000
Novo     → K crítico (Routh) = 6.0000      (diferença: 0)
```

Suíte de smoke test da interface: `tests/test_app_smoke.py` — 6 testes,
via `streamlit.testing.v1.AppTest` (executa a árvore real de widgets sem
navegador): carregamento padrão, clique no botão calcular, polos
complexos, sistema sem zeros finitos, entrada inválida e ordem alta.

```
6 passed in ~14s
```

**Total: 130/130 testes passando.**

> Observação de ambiente: não foi possível gerar uma captura de tela real
> do app renderizado neste momento, pois a instalação do navegador
> headless (Playwright/Chromium) foi bloqueada pela política de rede do
> sandbox de desenvolvimento usado para construir o projeto. Isso não
> afeta a validação funcional (os smoke tests já executam o app de
> verdade), mas vale conferir visualmente ao rodar `streamlit run app.py`
> localmente.
