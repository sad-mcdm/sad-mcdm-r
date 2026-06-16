# sadMCDM: Multi-Criteria Decision Making Solvers

Biblioteca R unificada de resolvedores e simulação de Monte Carlo para Apoio à Decisão Multicritério (MCDM - Multi-Criteria Decision Making). Esta biblioteca foi criada em R como o equivalente exato do pacote Python `sad-mcdm-lib`, fornecendo consistência matemática e suporte para o ecossistema SAD MCDM.

## Instalação

Para instalar a biblioteca localmente a partir do diretório local:

```R
# Requer o pacote devtools
install.packages("devtools")

# Instalação local
devtools::install_local("h:/Meu Drive/Pedro - CDSID/MCDM/sad-mcdm-r", force = TRUE)
```

Ou diretamente pelo GitHub (após subir o repositório):

```R
# Instalação pelo GitHub
devtools::install_github("pedrogouveia001/sad-mcdm-r")
```

## Dependências

O pacote utiliza pacotes padrão de R e requer:
- `lpSolve` (para os resolvedores baseados em programação linear: BWM, BWT, MACBETH, PROMETHEE V)
- `testthat` (opcional, para executar a suíte de testes unitários)

---

## Estrutura de Uso

Todas as funções são exportadas pelo namespace do pacote. Para utilizá-las, carregue o pacote no R:

```R
library(sadMCDM)
```

---

## Documentação das Funções e Resolvedores

Os argumentos em R seguem uma padronização limpa:
- `matrix_data`: Uma lista de listas contendo os valores físicos de consequências (por exemplo, `list(list(10, 100), list(20, 50), list(30, 0))`).
- `criteria_types`: Um vetor de strings contendo `"benefit"` (benefício) ou `"cost"` (custo).
- `preference_data`: Uma lista nomeada com parâmetros de preferência (pesos, limiares, matrizes de julgamento).
- `criteria_ids`: Um vetor contendo os IDs identificadores dos critérios (e.g. `c(1, 2)`).
- `alternatives_ids`: Um vetor contendo os IDs identificadores das alternativas (e.g. `c(1, 2, 3)`).

### 1. `solve_ahp`
Soluciona problemas de decisão usando o Processo Hierárquico Analítico (AHP) através do método da Média Geométrica das Linhas (RGM) e calcula a Razão de Consistência (CR).

* **Assinatura**:
  ```R
  solve_ahp(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `criteria_matrix`: Matriz de comparação de critérios $n \times n$.
  - `alternatives_matrices` (opcional): Lista contendo matrizes de comparações de alternativas para cada critério. Se não fornecida, as prioridades das alternativas serão derivadas automaticamente pela razão física das consequências.
* **Retorno**: Lista contendo:
  - `weights`: Pesos finais calculados para os critérios.
  - `normalized_matrix`: Matriz normalizada de prioridades locais das alternativas.
  - `global_scores`: Pontuações globais das alternativas.
  - `ranks`: Ordenação (ranks) das alternativas.
  - `criteria_cr`: Razão de Consistência (CR) da matriz de critérios.
  - `alternatives_cr`: Lista com a Razão de Consistência (CR) das alternativas por critério.

---

### 2. `solve_bwm`
Soluciona problemas de decisão usando o Método Best-Worst (BWM) através da resolução de um modelo de programação linear minimax.

* **Assinatura**:
  ```R
  solve_bwm(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `best_idx`: Índice (1-based) do melhor critério.
  - `worst_idx`: Índice (1-based) do pior critério.
  - `best_to_others`: Vetor de tamanho $n$ contendo a comparação do melhor critério com os outros.
  - `others_to_worst`: Vetor de tamanho $n$ contendo a comparação dos outros critérios com o pior.
  - `alternatives_matrices` (opcional): Idêntico ao do AHP para comparar alternativas.
* **Retorno**: Lista contendo `weights`, `normalized_matrix`, `global_scores`, `ranks`, `criteria_success`, `consistency_xi` (o valor do parâmetro $\xi^*$ de consistência) e `alternatives_cr`.

---

### 3. `solve_bwt`
Soluciona problemas de decisão utilizando o Best-Worst Tradeoff (BWT) com bisseção de funções de utilidade local por programação linear.

* **Assinatura**:
  ```R
  solve_bwt(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `best_idx`, `worst_idx`, `best_to_others`, `others_to_worst`: Mesmos julgamentos ordinais e cardinais inter-critérios do BWM.
  - `bisection_midpoints`: Lista nomeada associando cada critério ID com o ponto físico médio no qual a utilidade local é avaliada em exatamente 0.5.
* **Retorno**: Lista contendo `weights`, `normalized_matrix` (valores normalizados por interpolação de bisseção por partes), `global_scores` e `ranks`.

---

### 4. `solve_electre`
Soluciona problemas de decisão usando a família de métodos de sobreclassificação ELECTRE (suporta ELECTRE I, II, III, IV e TRI).

* **Assinatura**:
  ```R
  solve_electre(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `electre_version`: String indicando a versão (e.g. `"I"`, `"II"`, `"III"`, `"IV"`, `"TRI"`).
  - `weights`: Lista nomeada com os pesos numéricos de importância para cada critério.
  - `concordance_threshold` / `discordance_threshold`: Limiares para cálculo do grafo de sobreclassificação (versões I / IS).
  - `thresholds`: Lista contendo limiares de indiferença `q`, preferência `p` e veto `v` para cada critério ID.
* **Retorno**: Lista com `global_scores` (se aplicável), `ranks`, e um objeto `extra` com as relações de sobreclassificação (e.g., `kernel` ou classes preditas no ELECTRE TRI).

---

### 5. `solve_macbeth`
Soluciona problemas usando o método MACBETH (Measuring Attractiveness by a Categorical Based Evaluation Technique) através de programação linear baseada em julgamentos qualitativos semânticos.

* **Assinatura**:
  ```R
  solve_macbeth(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `criteria_matrix`: Matriz semântica qualitativa de atratividade para critérios.
  - `levels_matrices`: Lista nomeada de matrizes semânticas qualitativas de atratividade de 5 níveis de desempenho para interpolação de escalas locais de utilidade.
* **Retorno**: Lista contendo `weights`, `normalized_matrix`, `global_scores`, `ranks` e flags de sucesso (`criteria_success` e `levels_success`).

---

### 6. `solve_promethee`
Soluciona problemas de decisão usando os métodos PROMETHEE I (ordenação parcial), PROMETHEE II (ordenação completa baseada em fluxos líquidos) e PROMETHEE V (otimização de portfólio restrita).

* **Assinatura**:
  ```R
  solve_promethee(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `promethee_version`: Versão (`"I"`, `"II"` ou `"V"`).
  - `weights`: Lista nomeada contendo os pesos brutos dos critérios.
  - `thresholds`: Parâmetros de indiferença `q` e preferência `p` para cada critério.
  - No caso da versão `"V"`, deve incluir dados de custo por alternativa e restrição orçamentária máxima.
* **Retorno**: Lista com `phi_plus` (fluxo de saída positivo), `phi_minus` (fluxo de entrada negativo), `global_scores` (fluxo líquido líquido), `ranks` e dados extras do portfólio ótimo.

---

### 7. `solve_smarts_smarter`
Soluciona problemas usando as técnicas SMARTS (com pesos swing fornecidos diretamente) e SMARTER (calculando pesos centróides ROC a partir de uma ordenação ordinal de critérios).

* **Assinatura**:
  ```R
  solve_smarts_smarter(matrix_data, criteria_types, preference_data, criteria_ids, method)
  ```
* **Argumentos**:
  - `preference_data`: Lista contendo `weights` (pesos swing informados pelo decisor para o SMARTS) ou `ranks` (vetor com os IDs dos critérios ordenados do mais importante para o menos importante para o SMARTER).
  - `method`: String contendo `"smarts"` ou `"smarter"`.
* **Retorno**: Lista contendo `weights` (pesos calibrados ou ROC), `normalized_matrix`, `global_scores` e `ranks`.

---

### 8. `solve_topsis`
Soluciona problemas usando o método TOPSIS baseado na proximidade geométrica das alternativas em relação à solução ideal positiva (ideal) e negativa (anti-ideal).

* **Assinatura**:
  ```R
  solve_topsis(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `weights`: Lista nomeada com os pesos dos critérios (associando o critério ID ao peso, e.g. `list("1" = 40, "2" = 60)`).
* **Retorno**: Lista contendo `weights` normalizados, `normalized_matrix`, `global_scores` (coeficiente de proximidade relativa $C_i$), `ranks`, `distance_positive` (distância para ideal positivo) e `distance_negative` (distância para ideal negativo).

---

### 9. `solve_vikor`
Soluciona problemas usando o método de compromisso VIKOR, integrando a utilidade de grupo máxima e o arrependimento individual.

* **Assinatura**:
  ```R
  solve_vikor(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)
  ```
* **Argumentos em `preference_data`**:
  - `weights`: Lista nomeada com os pesos dos critérios.
  - `v`: Peso da estratégia de decisão da maioria (geralmente `0.5`).
* **Retorno**: Lista contendo `weights`, `normalized_matrix`, `global_scores` (compromisso global $Q_i$), `ranks`, `S` (distância da utilidade de grupo), `R` (arrependimento individual máximo) e `v`.

---

### 10. `run_monte_carlo` (Análise de Sensibilidade Estocástica)
Perturba a matriz física de consequências gerando ruídos uniformes dentro de uma faixa percentual por critério para avaliar a estabilidade do ranking final sob incerteza.

* **Assinatura**:
  ```R
  run_monte_carlo(matrix_data, criteria_types, weights, variations_pct, num_simulations, method, preference_data, criteria_ids)
  ```
* **Argumentos**:
  - `weights`: Vetor numérico contendo os pesos normalizados (soma = 1.0) para os critérios.
  - `variations_pct`: Vetor numérico com o percentual de variação tolerado em relação à amplitude de cada critério.
  - `num_simulations`: Quantidade de simulações (e.g. `100` ou `10000`).
  - `method`: String com o método que será perturbado (`"topsis"`, `"vikor"`, `"electre"`, `"promethee"`, `"bwt"`, `"macbeth"`, ou `"linear"`).
  - `preference_data` e `criteria_ids`: Parâmetros de preferência do resolvedor a ser executado.
* **Retorno**: Lista contendo:
  - `first_place_probabilities`: Vetor com a probabilidade de cada alternativa ficar em primeiro lugar.
  - `rank_probabilities`: Lista onde cada elemento é a probabilidade da alternativa $i$ ficar em cada uma das posições de rank de $1$ a $m$.
  - `average_scores`: Vetor com a pontuação média de cada alternativa.
  - `average_ranks`: Vetor com o ranking médio obtido por cada alternativa.
  - `deltas`: Amplitude máxima da variação aplicada a cada critério.

---

## Exemplo Rápido de Uso (TOPSIS)

Abaixo, um script simples exemplificando como chamar o resolvedor TOPSIS no R:

```R
library(sadMCDM)

# 3 alternativas, 2 critérios
matrix_data <- list(
  list(10.0, 100.0), # Alternativa 1
  list(20.0, 50.0),  # Alternativa 2
  list(30.0, 0.0)    # Alternativa 3
)

criteria_types <- c("benefit", "cost") # Critério 1 é benefício, Critério 2 é custo

preference_data <- list(
  weights = list("1" = 50.0, "2" = 50.0) # Pesos iguais para critério ID 1 e ID 2
)

criteria_ids <- c(1, 2)
alternatives_ids <- c(1, 2, 3)

# Executa o TOPSIS
resultado <- solve_topsis(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids)

# Exibe resultados
print("Pontuações (Closeness):")
print(resultado$global_scores)

print("Rankings Finais:")
print(resultado$ranks)
```
