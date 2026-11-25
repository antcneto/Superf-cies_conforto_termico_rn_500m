# Conforto térmico no Rio Grande do Norte: código de replicação

Este repositório reúne o script em R que reproduz as análises, as tabelas, as figuras e as superfícies de 500 m do artigo:

> CAMPOS NETO, A. Conforto térmico no Rio Grande do Norte: limites da interpolação de índices em rede escassa. **E&S Engineering and Science**, 2026. No prelo.

As superfícies interpoladas, a documentação metodológica e os guias de uso estão disponíveis em acesso aberto no Zenodo: <https://zenodo.org/records/23482528>.

## Sobre o estudo

O estudo avalia o que uma rede escassa de estações automáticas permite afirmar sobre o conforto térmico no Rio Grande do Norte e qual é o papel do vento na diferença entre o Índice de Calor (IC), que ignora o vento, e a Temperatura Efetiva Normal (NET), que o considera. Com dados horários de 2025 das redes da EMPARN e do INMET, quinze métodos de interpolação foram comparados por validação cruzada às 00h, 09h e 15h (hora local), tendo a reanálise ERA5-Land como referência. A diferença entre os índices foi decomposta em uma parcela do vento e outra da temperatura e da umidade, e as superfícies das 15h foram agregadas para os 167 municípios do estado. Às 15h, as regressões com altitude e reanálise tiveram os menores erros, e 31% dos 7,7 °C de diferença média entre os índices decorrem do vento.

## Conteúdo

| Arquivo | Descrição |
|---|---|
| `replicacao_conforto_termico_RN_2025.R` | Script único. O cabeçalho contém instruções de uso, glossário, lista dos dados de entrada e descrição das saídas. |

## O que o script reproduz

| Etapa | Conteúdo | No artigo |
|---|---|---|
| A | Leitura das observações, verificação da referência temporal de cada rede e controle de qualidade | Seção 2.2 |
| B | IC e NET por observação, redução do vento do INMET a 2 m e critério de completude (180 registros válidos) | Seções 2.2 e 2.3, Tabela 1, Figura 1 |
| C | Grade de 500 m (800 × 480 células, UTM zona 24S) e valores das covariáveis nas estações | Seção 2.4 |
| D | Validação cruzada *leave-one-out* dos métodos M1 a M15 e das referências R1 e R2 nos três horários | Seções 2.5 a 2.7, Quadro 1, Tabelas 2 e 3, Figura 2 |
| E | Intervalos de confiança do RMSE por *bootstrap* (5.000 réplicas) e testes pareados contra R2 | Seção 2.7, Tabela S2 |
| F | Diagnóstico da krigagem: variogramas por rodada, variograma fixo e REML | Seções 2.7 e 3.2, Tabela S1 |
| G | Decomposição da diferença entre IC e NET, classes de vento, ciclo diário e análises de sensibilidade | Seções 3.5 e 3.6, Tabela 4, Figuras 3 e 4 |
| H | NET recomposta a partir de T, UR e V interpolados | Seções 2.8 e 3.7, Tabela 5 |
| I | Validação cruzada sem o critério de completude (22 estações) | Seção 2.7 |
| J | Superfícies de 500 m, distância à estação mais próxima e zonas de confiança | Seções 2.9 e 3.8, Figura 5 |
| K | Médias municipais, população e classificação pelas escalas do NWS e de Fortaleza | Seções 2.9 e 3.8, Tabela 6, Figuras 6 e 7 |
| L | Comparação entre as duas estações de Natal | Seção 3.9 |
| M e N | Tabelas e figuras formatadas | Tabelas 1 a 6, S1 e S2, Figuras 1 a 7 |
| O | Conferência automática de 96 valores publicados no artigo | |

### Métodos avaliados

| ID | Método | Código |
|---|---|---|
| M1 | Inverso do quadrado da distância | `IDW` |
| M2 | Krigagem ordinária | `OK` |
| M3 | Krigagem com deriva: altitude | `KDE_alt` |
| M4 | Krigagem com deriva: altitude e NDVI | `KDE_alt_ndvi` |
| M5 | Regressão: coordenadas | `LM_coords` |
| M6 | Regressão: coordenadas e altitude | `LM_coords_alt` |
| M7 | Regressão: coordenadas, altitude e zona climática local | `LM_coords_alt_lcz` |
| M8 | GAM: coordenadas | `GAM_coords` |
| M9 | GAM: coordenadas e altitude | `GAM_coords_alt` |
| M10 | Regressão: altitude e reanálise | `LM_alt_era5` |
| M11 | Regressão: coordenadas, altitude e reanálise | `LM_coords_alt_era5` |
| M12 | Krigagem com deriva: altitude e reanálise | `KDE_alt_era5` |
| M13 | Floresta aleatória: coordenadas, altitude e reanálise | `RF_era5` |
| M14 | Regressão: altitude, reanálise e distância à costa | `LM_alt_era5_dist` |
| M15 | Krigagem com deriva: altitude, reanálise e distância à costa | `KDE_alt_era5_dist` |
| R1 | ERA5-Land no ponto da estação | `ERA5_bruto` |
| R2 | ERA5-Land corrigido pelo viés médio das demais estações | `ERA5_corrigido` |

## Requisitos

R 4.3 ou superior, com os pacotes:

```r
install.packages(c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2",
                   "sf", "sp", "terra", "gstat", "automap", "mgcv", "ranger",
                   "geobr", "hexbin", "patchwork", "ggrepel"))
```

Os quatro últimos são opcionais: o geobr baixa os limites do estado e dos municípios (dispensável quando eles são informados em arquivos `.gpkg` ou `.shp`), e hexbin, patchwork e ggrepel aprimoram as figuras. O pacote LCZ4r (Anjos et al. 2025) só é necessário para baixar o mapa de zonas climáticas locais.

Versões usadas no artigo:

* validação cruzada dos quinze métodos: R 4.3.3, terra 1.7.65, sf 1.0.15, sp 2.1.2, gstat 2.1.1, automap 1.1.20, mgcv 1.9.1 e ranger 0.16.0;
* análises complementares: R 4.6.1, terra 1.9.46, sf 1.1.2, sp 2.2.3, gstat 2.1.6, automap 1.1.20, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, purrr 1.2.2, ggplot2 4.0.3 e geobr 2.0.1.

Cada execução grava as versões efetivamente usadas em `15_conferencia/sessao_R.txt`.

## Como usar

1. Organize os dados de entrada conforme a seção seguinte.
2. Informe a pasta do projeto, editando `DIR_PROJETO` na seção 5 do script ou definindo a variável de ambiente `DIR_PROJETO_RN` antes da execução. Os demais caminhos e nomes de colunas também ficam na seção 5, único bloco a editar.
3. Execute o script inteiro numa sessão nova do R (use barras normais `/` também no Windows):

```r
Sys.setenv(DIR_PROJETO_RN = "C:/caminho/para/RN")
source("replicacao_conforto_termico_RN_2025.R", encoding = "UTF-8")
```

4. Ao final, abra `replicacao_artigo_2025/15_conferencia/conferencia_valores_publicados.csv`, que compara 96 valores publicados no artigo com os obtidos e classifica cada um como `confere`, `difere` ou `não calculado`.

A execução completa leva de 15 a 40 minutos, a maior parte nas validações cruzadas. As etapas demoradas gravam resultados intermediários em arquivos `.rds`, reaproveitados nas execuções seguintes enquanto as entradas não mudarem (`REAPROVEITAR = TRUE`).

## Dados de entrada

| Dado | Fonte | Como obter |
|---|---|---|
| Observações horárias de 2025 das 8 estações automáticas do INMET | INMET | <https://portal.inmet.gov.br/dadoshistoricos> (arquivo `2025.zip`); com `BAIXAR_INMET = TRUE`, o script baixa e descompacta o arquivo |
| Observações horárias de 2025 das 17 plataformas de coleta de dados da EMPARN | EMPARN | Não acompanham este repositório; podem ser solicitadas à EMPARN |
| Altitude (SRTM), NDVI (MOD13Q1 v061) e ERA5-Land horário | NASA e ECMWF, via Google Earth Engine | Código gravado pelo próprio script em `16_apendice/gee_covariaveis_RN_2025.js` |
| Zonas climáticas locais | Demuzere et al. (2022) | Pacote LCZ4r (`BAIXAR_LCZ_COM_LCZ4R = TRUE`) |
| Limites do estado e dos 167 municípios | IBGE | Pacote geobr (requer internet) ou arquivos `.gpkg` ou `.shp` |
| População residente no Censo 2022 | IBGE | SIDRA, tabela 4709: CSV com o código de 7 dígitos e a população de cada município |

As observações podem ser fornecidas de duas formas, escolhidas em `MODO_OBS`:

* `"padronizado"`: um único CSV, com uma linha por estação e hora e as colunas `estacao`, `rede` (`EMPARN` ou `INMET`), `lon`, `lat`, `datahora` (UTC, por exemplo `2025-01-01T03:00:00Z`), `T` (°C), `UR` (%) e `V` (m/s);
* `"brutos"`: os CSV oficiais do INMET e o arquivo da EMPARN. Nesse modo, o carimbo de tempo da EMPARN é tratado como hora local, conforme a verificação descrita na Seção 2.2 do artigo, e o script grava o arquivo padronizado.

Os arquivos do ERA5-Land devem conter as bandas `T_era5`, `UR_era5`, `V2m_era5`, `IC_era5` e `NET_era5_v2m`, geradas pelo código do Earth Engine.

Com os caminhos padrão, a pasta do projeto segue esta organização:

```
RN/                                                      pasta do projeto (DIR_PROJETO)
  artigo_vento/saidas/obs_horarias_2025.csv              observações padronizadas
  covariaveis/RN_srtm.tif                                altitude
  covariaveis/RN_ndvi_2025.tif                           NDVI médio de 2025
  covariaveis/lcz_rn.tif                                 zonas climáticas locais
  covariaveis/RN_era5land_h03_2025.tif                   ERA5-Land às 03 UTC (00h local)
  covariaveis/RN_era5land_h12_2025.tif                   ERA5-Land às 12 UTC (09h local)
  covariaveis/RN_era5land_h18_2025.tif                   ERA5-Land às 18 UTC (15h local)
  saidas/tabelas/municipios_classificados_populacao.csv  população por município
```

Arquivos opcionais:

```
  estacoes/inmet_2025/                                   CSV oficiais do INMET (modo "brutos")
  estacoes/emparn_rn_2025_unificado.csv                  dados da EMPARN (modo "brutos")
  estacoes/catalogo_estacoes.csv                         coordenadas das 25 estações (Figura 1)
  saidas/rasters/covariaveis_500m_hXX_2025.tif           pilha original de covariáveis
  artigo_vento/zenodo/RN_IC_15h_2025.tif                 grade de referência publicada no Zenodo
  artigo_vento/saidas/                                   saídas originais, para comparação
```

## Saídas

Os resultados são gravados em `replicacao_artigo_2025/`:

| Pasta | Conteúdo |
|---|---|
| `01_observacoes` | Observações padronizadas, controle de qualidade e referência temporal |
| `02_estacoes` | Médias por estação e horário e situação no critério de completude |
| `03_covariaveis` | Grade, máscaras, covariáveis e valores nas estações |
| `04_validacao` | Predições e métricas da validação cruzada e métodos vencedores |
| `05_incerteza` | Intervalos de *bootstrap* e testes pareados |
| `06_krigagem` | Parâmetros dos variogramas e variantes da krigagem |
| `07_decomposicao` | Decomposição, classes de vento, ciclo diário e sensibilidade |
| `08_recomposta` | NET recomposta |
| `09_completude` | Validação cruzada com 22 estações |
| `10_superficies` | GeoTIFFs de IC, NET, T, UR e V, distância às estações e zonas de confiança |
| `11_municipios` | Médias municipais, população e classes |
| `12_natal` | Comparação das estações de Natal |
| `13_tabelas` | Tabelas do artigo (CSV com ponto e vírgula e vírgula decimal) |
| `14_figuras` | Figuras do artigo (PNG, 300 dpi) |
| `15_conferencia` | Conferência com os valores publicados e informações da sessão do R |
| `16_apendice` | Código do Google Earth Engine |

## Reprodutibilidade

* Os índices, as médias por estação, a decomposição, as análises de sensibilidade e a comparação de Natal dependem apenas das observações e reproduzem o artigo com o mesmo arquivo de entrada.
* A validação cruzada e as superfícies dependem também das covariáveis. Para reproduzir os valores publicados, use a pilha original de covariáveis (`USAR_PILHA_ORIGINAL = TRUE`). Com `FALSE`, as covariáveis são refeitas a partir dos GeoTIFFs; esse preparo foi reescrito a partir da descrição dos Métodos, e as diferenças em relação aos valores originais ficam registradas em `03_covariaveis`.
* A floresta aleatória (M13) e as krigagens são sensíveis a diferenças mínimas de arredondamento nas entradas e podem variar na segunda ou na terceira casa decimal; os métodos vencedores, que são regressões lineares, não variam.
* O *bootstrap* usa a mesma semente do artigo (2025); os intervalos podem variar na segunda casa decimal se a ordem das estações for diferente.

## Como citar

> CAMPOS NETO, A. Conforto térmico no Rio Grande do Norte: limites da interpolação de índices em rede escassa. **E&S Engineering and Science**, 2026. No prelo.

Superfícies e documentação: <https://zenodo.org/records/23482528>.

## Licença

Código distribuído sob a licença MIT. Os dados de terceiros (INMET, EMPARN, IBGE, ECMWF e NASA) seguem as condições de uso de cada fonte.

## Autor

Antônio Campos Neto, Licenciado em Geografia, Universidade Federal do Rio Grande do Norte (UFRN).

Dúvidas e problemas podem ser registrados na aba *Issues* deste repositório.

## English summary

This repository contains a single R script that reproduces the analyses, tables, figures and 500 m surfaces of the article "Thermal comfort in Rio Grande do Norte, Brazil: limits of index interpolation over a sparse station network" (E&S Engineering and Science, 2026, in press). The study compares the Heat Index (HI) and the Normal Effective Temperature (NET) using 2025 hourly data from the EMPARN and INMET automatic station networks, evaluates fifteen interpolation methods by leave-one-out cross-validation against the ERA5-Land reanalysis, splits the difference between the two indices into a wind component and a temperature and humidity component, and aggregates the 3 p.m. surfaces to the 167 municipalities of the state. Code comments and outputs are in Portuguese; the script header lists every input dataset required.
