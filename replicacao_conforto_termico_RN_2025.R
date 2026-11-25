# =============================================================================
#  REPLICAÇÃO INTEGRAL DO ARTIGO
#  "Conforto térmico no Rio Grande do Norte: limites da interpolação de índices
#   em rede escassa" (Campos Neto, ES Engineering and Science)
#
#  Arquivo único em R que reproduz, a partir dos dados de entrada, as análises,
#  as tabelas, as figuras e as superfícies do artigo.
#  Versão 1.0, novembro de 2026.
# =============================================================================
#
#  SUMÁRIO
#    0. Como usar
#    1. Glossário
#    2. Dados necessários de antemão
#    3. O que o script produz
#    4. Pacotes
#    5. Configuração (único bloco a editar)
#    Etapa A  Observações: leitura, referência temporal e controle de qualidade
#    Etapa B  Índices por observação e médias por estação e horário
#    Etapa C  Limites, grade de 500 m e covariáveis
#    Etapa D  Validação cruzada (M1 a M15, R1 e R2) nos três horários
#    Etapa E  Incerteza do RMSE: bootstrap e testes pareados
#    Etapa F  Diagnóstico da krigagem (Tabela S1)
#    Etapa G  Decomposição IC - NET, classes de vento, ciclo diário e sensibilidade
#    Etapa H  NET recomposta a partir das variáveis interpoladas
#    Etapa I  Efeito do critério de completude (22 estações)
#    Etapa J  Superfícies de 500 m, distância às estações e zonas de confiança
#    Etapa K  Leitura municipal, população e escalas de classificação
#    Etapa L  Estações de Natal (INMET e EMPARN)
#    Etapa M  Tabelas do artigo
#    Etapa N  Figuras do artigo
#    Etapa O  Conferência com os valores publicados
#    Apêndice Código do Google Earth Engine que gera as covariáveis
#
# -----------------------------------------------------------------------------
#  0. COMO USAR
# -----------------------------------------------------------------------------
#  1. Instale uma vez os pacotes listados na seção 4.
#  2. Organize os arquivos de entrada conforme a seção 2. Na pasta do projeto do
#     autor, os caminhos padrão da seção 5 já apontam para os arquivos existentes.
#  3. Edite SOMENTE o bloco da seção 5 (pastas, nomes de arquivos e de colunas).
#  4. Rode o arquivo inteiro, de preferência numa sessão nova do R:
#       source("replicacao_conforto_termico_RN_2025.R", encoding = "UTF-8")
#     ou, no RStudio, Ctrl + Shift + Enter (Source).
#  5. Ao final, abra 15_conferencia/conferencia_valores_publicados.csv, que
#     compara os números obtidos com os publicados no artigo.
#
#  Tempo estimado: de 15 a 40 minutos, quase todo nas validações cruzadas
#  (Etapas D e I). As etapas demoradas gravam seus resultados em arquivos .rds;
#  com REAPROVEITAR = TRUE, uma nova execução os reutiliza. Para refazer tudo,
#  use REAPROVEITAR = FALSE ou apague a pasta de saída.
#
#  Reprodutibilidade:
#  - Os índices, as médias por estação, a decomposição, as análises de
#    sensibilidade e a comparação de Natal dependem só das observações e
#    reproduzem o artigo com o mesmo arquivo de entrada.
#  - A validação cruzada e as superfícies dependem também das covariáveis nas
#    estações. O preparo das covariáveis (reamostragem para 500 m, distância à
#    costa e valores nas estações) foi reescrito a partir da descrição dos
#    Métodos, porque o código que gerou a pilha original não foi preservado.
#    Para reproduzir exatamente os números publicados, use a pilha original
#    (USAR_PILHA_ORIGINAL = TRUE); para refazer tudo a partir dos dados brutos,
#    use FALSE e confira as diferenças em 03_covariaveis/comparacao_*.csv.
#  - A floresta aleatória (M13) e as krigagens são sensíveis a diferenças de
#    arredondamento nas entradas e podem variar na segunda ou terceira casa
#    decimal; os métodos vencedores (regressões lineares) não.
#  - Os intervalos de bootstrap dependem da sequência de números aleatórios e
#    da ordem das estações; podem variar na segunda casa decimal se essa ordem
#    diferir da original. A semente é a mesma do artigo (2025).
#
# -----------------------------------------------------------------------------
#  1. GLOSSÁRIO
# -----------------------------------------------------------------------------
#  IC        Índice de Calor (Rothfusz 1990, com os ajustes do NWS), em °C.
#  NET       Temperatura Efetiva Normal (Li e Chan 2000), em °C; usa T, UR e V.
#  NET0      NET calculada com vento nulo (ar parado), usada na decomposição.
#  T, UR, V  Temperatura do ar (°C), umidade relativa (%) e velocidade do vento (m/s).
#  V2, V10   Vento a 2 m e a 10 m. A EMPARN mede a 2 m e o INMET, a 10 m.
#  z0        Comprimento de rugosidade da lei logarítmica do vento (0,03 m).
#            Fator de redução de 10 m para 2 m: ln(2/z0) / ln(10/z0) = 0,723.
#  Decomposição  IC - NET = (IC - NET0) + (NET0 - NET): a primeira parcela é a de
#            temperatura e umidade; a segunda, a do vento.
#  EMPARN    Empresa de Pesquisa Agropecuária do RN (17 plataformas de coleta de dados, PCD).
#  INMET     Instituto Nacional de Meteorologia (8 estações automáticas no RN).
#  Hora local  UTC menos 3 horas. Horários analisados: 00h, 09h e 15h locais
#            (03, 12 e 18 UTC).
#  Critério de completude  Mínimo de 180 registros válidos no ano, no horário.
#  Registro válido  T, UR e V válidas simultaneamente após o controle de qualidade
#            (5 a 50 °C, 5 a 100% e 0 a 30 m/s).
#  LOOCV     Validação cruzada leave-one-out: ajusta com todas as estações menos uma
#            e prediz a excluída, para cada estação.
#  RMSE      Raiz do erro quadrático médio da validação cruzada.
#  R²        1 - SQE/SQT; negativo quando a predição é pior que a média das estações.
#  M1 a M15  Os quinze métodos de interpolação do Quadro 1 (lista na Etapa D).
#  R1, R2    Referências da reanálise: ERA5-Land no ponto da estação, bruto (R1)
#            e corrigido pelo viés médio das demais estações (R2).
#  IDW       Inverso do quadrado da distância (potência 2).
#  Krigagem  Interpolação geoestatística baseada no variograma. Ordinária (OK) ou
#            com deriva externa (KDE), quando há covariáveis na tendência.
#  Variograma  Função que descreve como a semelhança entre pontos cai com a
#            distância. Parâmetros: efeito pepita, patamar parcial e alcance.
#  REML      Máxima verossimilhança restrita (ajuste alternativo do variograma).
#  GAM       Modelo aditivo generalizado (mgcv), ajustado por REML.
#  RF        Floresta aleatória (ranger): 500 árvores, nó mínimo 2, mtry 2, semente 1.
#  ERA5-Land Reanálise do ECMWF (~9 km), extraída no Google Earth Engine.
#  NDVI      Índice de vegetação por diferença normalizada (MODIS MOD13Q1 v061).
#  LCZ       Zona climática local (mapa global de Demuzere et al. 2022, 100 m).
#  SRTM      Modelo digital de elevação da Shuttle Radar Topography Mission.
#  Grade     Grade de cerca de 500 m (800 x 480 células) em WGS 84 / UTM 24S.
#  Zonas de confiança  Faixas da distância de cada célula à estação mais próxima,
#            delimitadas pelos percentis 75 e 90 dessa distância.
#  Bootstrap Reamostragem das estações com reposição (5.000 réplicas) para o
#            intervalo de confiança de 95% do RMSE.
#  Wilcoxon  Teste pareado sobre os erros absolutos por estação (método contra R2).
#  NET recomposta  NET recalculada com T e UR preditos pelos métodos vencedores e
#            três origens do vento: rede, ERA5-Land a 2 m e vento nulo.
#  Escala NWS  Faixas do IC do National Weather Service: cuidado (27 a 32 °C),
#            cuidado extremo (32 a 41 °C), perigo (41 a 54 °C).
#  Escala de Fortaleza  Faixas de Lima Júnior et al. (2025) a partir de Petalas
#            (2015): confortável (29,3 a 32,1 °C), desconfortável (acima de 32,1 °C).
#
# -----------------------------------------------------------------------------
#  2. DADOS NECESSÁRIOS DE ANTEMÃO
# -----------------------------------------------------------------------------
#  A. OBSERVAÇÕES HORÁRIAS DAS ESTAÇÕES DE 2025 (obrigatório; escolha uma opção)
#
#   Opção 1, MODO_OBS = "padronizado": um único CSV, uma linha por estação e hora.
#     Colunas (os nomes podem ser ajustados em COLS_OBS):
#       estacao   nome da estação (texto)
#       rede      "EMPARN" ou "INMET"
#       lon, lat  coordenadas em graus decimais (WGS 84 ou SIRGAS 2000)
#       datahora  instante em UTC verdadeiro, p. ex. 2025-01-01T03:00:00Z
#       T, UR, V  temperatura (°C), umidade relativa (%) e vento (m/s)
#       altura_vento (opcional) altura do anemômetro em metros
#     O arquivo artigo_vento/saidas/obs_horarias_2025.csv já está nesse formato,
#     com o vento do INMET reduzido a 2 m (use V_JA_A_2M = TRUE). O arquivo que o
#     próprio script grava na Opção 2 (obs_horarias_padronizadas_2025.csv) também,
#     com o vento na altura nativa (use V_JA_A_2M = FALSE). Como os dados da EMPARN
#     não estão em acesso aberto, é este o arquivo a depositar no repositório.
#
#   Opção 2, MODO_OBS = "brutos":
#     (a) INMET: pasta com os CSV oficiais anuais das estações automáticas
#         (https://portal.inmet.gov.br/dadoshistoricos, arquivo 2025.zip). Com
#         BAIXAR_INMET = TRUE o script baixa e descompacta o arquivo sozinho.
#     (b) EMPARN: um CSV das PCDs com nome da estação, coordenadas, carimbo de data
#         e hora, T, UR, vento a 2 m e, de preferência, radiação global (usada na
#         verificação da referência temporal). Os nomes das colunas são informados
#         em COLS_EMPARN. O carimbo da EMPARN é tratado como hora local
#         (CARIMBO_EMPARN_EM_HORA_LOCAL = TRUE), conforme a Seção 2.2 do artigo.
#
#   Catálogo de estações (opcional): CSV com estacao, rede, lon, lat das 25
#     estações instaladas. Só é necessário para que a Figura 1 mostre as estações
#     sem nenhum registro válido (Apodi, Calcanhar e Lagoa Nova), caso elas não
#     constem do arquivo de observações.
#
#  B. COVARIÁVEIS (GeoTIFF em WGS 84, EPSG:4326). Todas, menos a LCZ, saem do
#     código do Google Earth Engine do apêndice (escrito em 16_apendice/*.js):
#       RN_srtm.tif              altitude (m), cerca de 90 m
#       RN_ndvi_2025.tif         NDVI médio de 2025 (MOD13Q1 v061), 250 m
#       lcz_rn.tif               LCZ (classes 1 a 17), 100 m, Demuzere et al. (2022),
#                                obtida com o pacote LCZ4r (BAIXAR_LCZ_COM_LCZ4R)
#       RN_era5land_h03_2025.tif médias anuais às 03 UTC (00h local)
#       RN_era5land_h12_2025.tif médias anuais às 12 UTC (09h local)
#       RN_era5land_h18_2025.tif médias anuais às 18 UTC (15h local)
#     Os três arquivos do ERA5-Land devem ter as bandas T_era5, UR_era5, V2m_era5,
#     IC_era5 e NET_era5_v2m, calculadas hora a hora e depois promediadas.
#
#  C. LIMITES: estado do RN, 167 municípios e linha de costa. Por padrão vêm do
#     pacote geobr (precisa de internet); também podem ser arquivos .gpkg ou .shp.
#     A linha de costa, se não for informada, é o limite do RN menos as divisas
#     com o Ceará e a Paraíba.
#
#  D. POPULAÇÃO: CSV com o código IBGE de 7 dígitos e a população residente do
#     Censo 2022 de cada município do RN (SIDRA, tabela 4709). Nomes das colunas
#     em COLS_POP.
#
#  E. OPCIONAIS:
#     Pilha original de covariáveis (saidas/rasters/covariaveis_500m_h03_2025.tif,
#                       h12 e h18), usada se USAR_PILHA_ORIGINAL = TRUE. Deve ter
#                       as bandas altitude, dist_costa (km), ndvi, lcz e as cinco
#                       bandas do ERA5-Land listadas acima.
#     ARQ_GRADE_MODELO  um GeoTIFF publicado no Zenodo (p. ex. RN_IC_15h_2025.tif).
#                       Se existir, a grade de 800 x 480 células é copiada dele, o
#                       que garante a mesma geometria das superfícies publicadas.
#     DIR_ORIGINAL      pasta com as saídas originais (artigo_vento/saidas). Se
#                       existir, as médias e as covariáveis das estações são
#                       comparadas com as originais.
#
# -----------------------------------------------------------------------------
#  3. O QUE O SCRIPT PRODUZ (dentro de DIR_SAIDA)
# -----------------------------------------------------------------------------
#  01_observacoes    observações padronizadas, controle de qualidade, referência temporal
#  02_estacoes       médias por estação e horário, situação no critério de completude
#  03_covariaveis    grade, máscaras, pilhas de covariáveis e valores nas estações
#  04_validacao      predições e métricas da LOOCV, vencedores (Tabelas 2 e 3)
#  05_incerteza      intervalos de bootstrap (Tabela S2) e testes pareados contra R2
#  06_krigagem       parâmetros dos variogramas por rodada (Tabela S1) e variantes
#  07_decomposicao   decomposição, classes de vento (Tabela 4), ciclo diário,
#                    sensibilidade, frequência de limiares, Jundiaí e derivadas
#  08_recomposta     NET recomposta (Tabela 5)
#  09_completude     validação com 22 estações e contraste com 16
#  10_superficies    GeoTIFFs de IC, NET, T, UR e V, distância e zonas de confiança
#  11_municipios     médias municipais, população, classes e deslocamentos
#  12_natal          comparação das duas estações de Natal
#  13_tabelas        Tabelas 1 a 6, S1 e S2 formatadas (CSV com ; e vírgula decimal)
#  14_figuras        Figuras 1 a 7 (PNG, 300 dpi), legíveis em tons de cinza
#  15_conferencia    comparação com os valores publicados e informações da sessão
#  16_apendice       código do Google Earth Engine (.js)
#
# -----------------------------------------------------------------------------
#  4. PACOTES
# -----------------------------------------------------------------------------
#  Obrigatórios: dplyr, tidyr, readr, purrr, tibble, ggplot2, sf, sp, terra,
#                gstat, automap, mgcv, ranger
#  Opcionais:    geobr (limites), hexbin (Figura 3a), patchwork (painéis),
#                ggrepel (rótulos da Figura 1), LCZ4r (download da LCZ)
#  Instalação:
#    install.packages(c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2",
#                       "sf", "sp", "terra", "gstat", "automap", "mgcv", "ranger",
#                       "geobr", "hexbin", "patchwork", "ggrepel"))
#  Versões do ambiente em que a validação cruzada do artigo foi calculada:
#    R 4.3.3; terra 1.7.65; sf 1.0.15; sp 2.1.2; gstat 2.1.1; automap 1.1.20;
#    mgcv 1.9.1; ranger 0.16.0. Análises complementares: R 4.6.1.
# =============================================================================


## ============================================================================
## 5. CONFIGURAÇÃO (único bloco a editar)
## ============================================================================
# Use barras normais "/" também no Windows.
DIR_PROJETO <- Sys.getenv("DIR_PROJETO_RN", "C:/Users/nirma/OneDrive/Documentos/Antonio/RN")
DIR_SAIDA   <- file.path(DIR_PROJETO, "replicacao_artigo_2025") # criada se não existir
# Os caminhos abaixo seguem a organização das pastas do projeto do autor. Em outra
# máquina, aponte cada um para o arquivo correspondente da seção 2.

# ---- A. Observações ---------------------------------------------------------
MODO_OBS <- "padronizado"            # "padronizado" (Opção 1) ou "brutos" (Opção 2)

# Opção 1
ARQ_OBS_PADRONIZADO <- file.path(DIR_PROJETO, "artigo_vento", "saidas", "obs_horarias_2025.csv")
COLS_OBS <- c(estacao = "estacao", rede = "rede", lon = "lon", lat = "lat",
              datahora = "datahora", T = "T", UR = "UR", V = "V")
V_JA_A_2M <- TRUE    # TRUE se o vento do INMET no arquivo já estiver reduzido a 2 m

# Opção 2
DIR_INMET_BRUTO  <- file.path(DIR_PROJETO, "estacoes", "inmet_2025")
BAIXAR_INMET     <- FALSE
ARQ_EMPARN_BRUTO <- file.path(DIR_PROJETO, "estacoes", "emparn_rn_2025_unificado.csv")
COLS_EMPARN <- list(estacao = "estacao", lon = "lon", lat = "lat",
                    datahora = "datahora",       # ou c("data", "hora") se separadas
                    T = "temperatura", UR = "umidade", V = "vento_2m_ms",
                    RAD = "radiacao")            # NA se não houver radiação
CARIMBO_EMPARN_EM_HORA_LOCAL <- TRUE

ARQ_CATALOGO <- file.path(DIR_PROJETO, "estacoes", "catalogo_estacoes.csv")  # opcional

# ---- B. Covariáveis -----------------------------------------------------------
DIR_COV  <- file.path(DIR_PROJETO, "covariaveis")
ARQ_SRTM <- file.path(DIR_COV, "RN_srtm.tif")
ARQ_NDVI <- file.path(DIR_COV, "RN_ndvi_2025.tif")
ARQ_LCZ  <- file.path(DIR_COV, "lcz_rn.tif")
ARQ_ERA5 <- c("0"  = file.path(DIR_COV, "RN_era5land_h03_2025.tif"),
              "9"  = file.path(DIR_COV, "RN_era5land_h12_2025.tif"),
              "15" = file.path(DIR_COV, "RN_era5land_h18_2025.tif"))
BAIXAR_LCZ_COM_LCZ4R <- FALSE        # TRUE: baixa a LCZ com LCZ4r se lcz_rn.tif não existir

# Pilha de covariáveis já reamostrada que alimentou o artigo (opcional). Com TRUE,
# as covariáveis não são recalculadas: são lidas destes arquivos (%02d = hora UTC).
USAR_PILHA_ORIGINAL <- FALSE
ARQ_PILHA_ORIGINAL  <- file.path(DIR_PROJETO, "saidas", "rasters", "covariaveis_500m_h%02d_2025.tif")

# ---- C e D. Limites e população -----------------------------------------------
ARQ_LIMITE_RN   <- NA                # NA = geobr; ou caminho de .gpkg/.shp
ARQ_MUNICIPIOS  <- NA                # NA = geobr; ou caminho de .gpkg/.shp
ARQ_LINHA_COSTA <- NA                # NA = derivada do limite do RN (geobr)
ANO_MALHA       <- 2020              # ano das malhas do geobr (confira as células na Etapa O)
ARQ_POPULACAO   <- file.path(DIR_PROJETO, "saidas", "tabelas", "municipios_classificados_populacao.csv")
COLS_POP        <- c(codigo = "code_muni", populacao = "populacao")   # outros nomes usuais são detectados

# ---- E. Opcionais ---------------------------------------------------------------
ARQ_GRADE_MODELO <- file.path(DIR_PROJETO, "artigo_vento", "zenodo", "RN_IC_15h_2025.tif")
DIR_ORIGINAL     <- file.path(DIR_PROJETO, "artigo_vento", "saidas")

# ---- Execução -------------------------------------------------------------------
REAPROVEITAR  <- TRUE                # reutiliza os .rds das etapas demoradas (se as entradas não mudaram)
GERAR_FIGURAS <- TRUE
## ======================= fim do bloco de configuração =======================


## ============================================================================
## Parâmetros fixados no artigo (não altere para reproduzir os resultados)
## ============================================================================
P <- list(
  ano = 2025, fuso = -3,                       # hora local = UTC - 3 h
  horas = c(0, 9, 15), hora_mapa = 15,
  min_registros = 180,
  qc = list(T = c(5, 50), UR = c(5, 100), V = c(0, 30)),
  z0 = 0.03, z0_sensibilidade = c(0.01, 0.10, 0.25), v_ref_sensibilidade = c(0.1, 0.2),
  crs = "EPSG:32724", ncol = 800, nrow = 480,
  raios_busca = c(1000, 3000, 6000, 12000),   # metros, para estações sem valor na célula
  reamostragem_continua = "bilinear", reamostragem_lcz = "mode",
  B = 5000, semente_bootstrap = 2025, semente_geral = 20261019,
  classes_vento = c(-Inf, 0.5, 1.5, 3, 4, 5, Inf),
  limiares_nws = c(27, 32, 41, 54), limiares_fortaleza = c(29.3, 32.1),
  limiares_frequencia = c(32, 41))

VARIAVEIS <- c("IC", "NET", "T", "UR", "V")
COL_ERA5  <- c(IC = "IC_era5", NET = "NET_era5_v2m", T = "T_era5",
               UR = "UR_era5", V = "V2m_era5")
UNIDADE   <- c(IC = "°C", NET = "°C", T = "°C", UR = "%", V = "m/s")

METODOS <- tibble::tribble(
  ~id,   ~codigo,              ~familia,            ~rotulo,
  "M1",  "IDW",                "Determinística",    "Inverso do quadrado da distância",
  "M2",  "OK",                 "Geoestatística",    "Krigagem ordinária",
  "M3",  "KDE_alt",            "Geoestatística",    "Krigagem com deriva: altitude",
  "M4",  "KDE_alt_ndvi",       "Geoestatística",    "Krigagem com deriva: altitude e NDVI",
  "M5",  "LM_coords",          "Regressão linear",  "Regressão: coordenadas",
  "M6",  "LM_coords_alt",      "Regressão linear",  "Regressão: coordenadas e altitude",
  "M7",  "LM_coords_alt_lcz",  "Regressão linear",  "Regressão: coordenadas, altitude e zona climática local",
  "M8",  "GAM_coords",         "GAM",               "GAM: coordenadas",
  "M9",  "GAM_coords_alt",     "GAM",               "GAM: coordenadas e altitude",
  "M10", "LM_alt_era5",        "Reanálise",         "Regressão: altitude e reanálise",
  "M11", "LM_coords_alt_era5", "Reanálise",         "Regressão: coordenadas, altitude e reanálise",
  "M12", "KDE_alt_era5",       "Reanálise",         "Krigagem com deriva: altitude e reanálise",
  "M13", "RF_era5",            "Reanálise",         "Floresta aleatória: coordenadas, altitude e reanálise",
  "M14", "LM_alt_era5_dist",   "Distância à costa", "Regressão: altitude, reanálise e distância à costa",
  "M15", "KDE_alt_era5_dist",  "Distância à costa", "Krigagem com deriva: altitude, reanálise e distância à costa")
REFERENCIAS <- tibble::tibble(id = c("R1", "R2"), codigo = c("ERA5_bruto", "ERA5_corrigido"),
                              familia = "Referência",
                              rotulo = c("ERA5-Land bruto", "ERA5-Land corrigido do viés"))
ROTULOS <- dplyr::bind_rows(METODOS, REFERENCIAS)

# Covariáveis exigidas por cada método (para descartar células ou estações incompletas)
NECESSARIAS <- list(
  IDW = c("x", "y"), OK = c("x", "y"), KDE_alt = c("x", "y", "altitude"),
  KDE_alt_ndvi = c("x", "y", "altitude", "ndvi"),
  LM_coords = c("x", "y"), LM_coords_alt = c("x", "y", "altitude"),
  LM_coords_alt_lcz = c("x", "y", "altitude", "lcz"),
  GAM_coords = c("x", "y"), GAM_coords_alt = c("x", "y", "altitude"),
  LM_alt_era5 = c("altitude", "era5"), LM_coords_alt_era5 = c("x", "y", "altitude", "era5"),
  KDE_alt_era5 = c("x", "y", "altitude", "era5"), RF_era5 = c("x", "y", "altitude", "era5"),
  LM_alt_era5_dist = c("altitude", "era5", "dist_costa"),
  KDE_alt_era5_dist = c("x", "y", "altitude", "era5", "dist_costa"))


## ============================================================================
## Pacotes e funções auxiliares
## ============================================================================
obrig <- c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2", "sf", "sp",
           "terra", "gstat", "automap", "mgcv", "ranger")
falta_pac <- obrig[!vapply(obrig, requireNamespace, logical(1), quietly = TRUE)]
if (length(falta_pac)) {
  stop("Instale antes de rodar: install.packages(c(",
       paste0('"', falta_pac, '"', collapse = ", "), "))", call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(tibble)
  library(ggplot2); library(sf); library(terra); library(mgcv)
})
options(dplyr.summarise.inform = FALSE, readr.show_progress = FALSE, warn = 1)
terra::terraOptions(progress = 0)
tem <- function(p) requireNamespace(p, quietly = TRUE)
set.seed(P$semente_geral)

msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
pasta <- function(nome) {
  d <- file.path(DIR_SAIDA, nome)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
gravar <- function(df, dir, nome) {
  readr::write_csv(df, file.path(dir, nome), na = "")
  invisible(df)
}
gravar_br <- function(df, dir, nome) {               # para abrir no Excel em português
  readr::write_excel_csv2(df, file.path(dir, nome), na = "")
  invisible(df)
}
# Guarda o resultado de uma etapa demorada; só o reutiliza se a assinatura das
# entradas (soma das colunas numéricas) for a mesma da execução anterior.
assinatura <- function(...) {
  tabs <- list(...)
  format(sum(vapply(tabs, function(d) {
    d <- as.data.frame(d)
    sum(vapply(d[vapply(d, is.numeric, logical(1))], function(x) sum(x, na.rm = TRUE), numeric(1))) + nrow(d)
  }, numeric(1))), digits = 15)
}
em_cache <- function(arq, expr, ass = NULL) {
  if (REAPROVEITAR && file.exists(arq)) {
    guardado <- readRDS(arq)
    if (is.list(guardado) && identical(names(guardado), c("valor", "assinatura")) &&
        (is.null(ass) || identical(guardado$assinatura, ass))) {
      msg("   reaproveitando ", basename(arq))
      return(guardado$valor)
    }
  }
  val <- expr
  saveRDS(list(valor = val, assinatura = ass), arq)
  val
}
existe <- function(f) length(f) == 1 && !is.na(f) && nzchar(f) && file.exists(f)
sem_acento <- function(x) {
  chartr("ÁÀÂÃÄÉÈÊËÍÌÎÏÓÒÔÕÖÚÙÛÜÇÑáàâãäéèêëíìîïóòôõöúùûüçñ",
         "AAAAAEEEEIIIIOOOOOUUUUCNaaaaaeeeeiiiiooooouuuucn", x)
}
chave <- function(x) toupper(trimws(gsub("\\s+", " ", sem_acento(as.character(x)))))
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
moda <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  t <- table(x)
  as.numeric(names(t)[which.max(t)])
}
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
num_br <- function(x, d = 1) {
  s <- formatC(x, format = "f", digits = d, big.mark = ".", decimal.mark = ",")
  s <- sub("^-", "\u2212", s)
  s[is.na(x)] <- ""
  s
}
fmt_r2 <- function(x) ifelse(is.na(x), "", ifelse(x < -10, "< \u221210", num_br(x, 3)))
fmt_rmse <- function(x, var) {
  var <- rep_len(var, length(x))
  vapply(seq_along(x), function(i) {
    v <- x[i]
    if (is.na(v)) return("")
    d <- if (v >= 100) 0 else if (v >= 10) 1 else if (var[i] == "UR") 2 else 3
    num_br(v, d)
  }, character(1))
}
ordenar_estacoes <- function(df) df[order(df$estacao), , drop = FALSE]
eh_jundiai <- function(x) grepl("JUNDIA", chave(x))

# Nomes para tabelas e figuras (chave sem acento -> rótulo)
NOMES <- tibble::tribble(
  ~chave,                            ~tabela,                         ~figura,
  "CAICO (A316)",                    "Caicó (A316)",                  "Caicó (A316)",
  "IPANGUACU (A372)",                "Ipanguaçu (A372)",              "Ipanguaçu (A372)",
  "MACAU (A317)",                    "Macau (A317)",                  "Macau (A317)",
  "MOSSORO (A318)",                  "Mossoró (A318)",                "Mossoró (A318)",
  "NATAL (A304)",                    "Natal (A304)",                  "Natal (A304)",
  "APODI (A340)",                    "Apodi (A340)",                  "Apodi (A340)",
  "CALCANHAR (A344)",                "Calcanhar (A344)",              "Calcanhar (A344)",
  "SANTA CRUZ (A367)",               "Santa Cruz (A367)",             "Santa Cruz (A367)",
  "PCD ACARI",                       "PCD Acari",                     "Acari",
  "PCD ALEXANDRIA",                  "PCD Alexandria",                "Alexandria",
  "PCD BOA SAUDE",                   "PCD Boa Saúde",                 "Boa Saúde",
  "PCD CANGUARETAMA",                "PCD Canguaretama",              "Canguaretama",
  "PCD CARAUBAS",                    "PCD Caraúbas",                  "Caraúbas",
  "PCD ESCOLA AGRICOLA DE JUNDIAI",  "PCD Esc. Agrícola de Jundiaí",  "Esc. Agr. de Jundiaí",
  "PCD GALINHOS",                    "PCD Galinhos",                  "Galinhos",
  "PCD NATAL",                       "PCD Natal",                     "Natal (EMPARN)",
  "PCD PUREZA",                      "PCD Pureza",                    "Pureza",
  "PCD SERRA DO MEL",                "PCD Serra do Mel",              "Serra do Mel",
  "PCD SAO GONCALO DO AMARANTE",     "PCD São Gonçalo do Amarante",   "São Gonçalo do Amarante",
  "PCD BARAUNA",                     "PCD Baraúna",                   "Baraúna",
  "PCD JOAO CAMARA",                 "PCD João Câmara",               "João Câmara",
  "PCD SANTANA DO MATOS",            "PCD Santana do Matos",          "Santana do Matos",
  "PCD MARTINS",                     "PCD Martins",                   "Martins",
  "PCD SAO MIGUEL",                  "PCD São Miguel",                "São Miguel",
  "PCD LAGOA NOVA",                  "PCD Lagoa Nova",                "Lagoa Nova")
rotulo_estacao <- function(x, tipo = c("tabela", "figura")) {
  tipo <- match.arg(tipo)
  r <- NOMES[[tipo]][match(chave(x), NOMES$chave)]
  ifelse(is.na(r), as.character(x), r)
}

# ---- Índices (iguais aos do pipeline do artigo) ----------------------------------
# IC: regressão de Rothfusz (1990) em °F, ajustes do NWS para UR < 13% e UR > 85%
# e forma simplificada quando a média entre ela e a temperatura fica abaixo de 80 °F.
ic_c <- function(T_c, RH) {
  Tf <- T_c * 9 / 5 + 32
  simples <- 0.5 * (Tf + 61 + (Tf - 68) * 1.2 + RH * 0.094)
  cheio <- -42.379 + 2.04901523 * Tf + 10.14333127 * RH - 0.22475541 * Tf * RH -
    0.00683783 * Tf^2 - 0.05481717 * RH^2 + 0.00122874 * Tf^2 * RH +
    0.00085282 * Tf * RH^2 - 0.00000199 * Tf^2 * RH^2
  aj1 <- ifelse(RH < 13 & Tf >= 80 & Tf <= 112,
                -((13 - RH) / 4) * sqrt(pmax(17 - abs(Tf - 95), 0) / 17), 0)
  aj2 <- ifelse(RH > 85 & Tf >= 80 & Tf <= 87, ((RH - 85) / 10) * ((87 - Tf) / 5), 0)
  hi <- ifelse((simples + Tf) / 2 < 80, simples, cheio + aj1 + aj2)
  (hi - 32) * 5 / 9
}
# NET: Li e Chan (2000); V em m/s a 2 m.
net_c <- function(T_c, RH, V) {
  37 - (37 - T_c) / (0.68 - 0.0014 * RH + 1 / (1.76 + 1.4 * V^0.75)) -
    0.29 * T_c * (1 - RH / 100)
}
# Redução do vento da altura h para 2 m pela lei logarítmica (Equação 3).
fator_altura <- function(h, z0 = P$z0) log(2 / z0) / log(h / z0)

# ---- Leitura de datas e horas ------------------------------------------------------
ler_datahora <- function(x) {
  if (inherits(x, "POSIXct")) { attr(x, "tzone") <- "UTC"; return(x) }
  if (inherits(x, "Date")) return(as.POSIXct(format(x), tz = "UTC"))
  s <- trimws(as.character(x))
  s <- sub("Z$", "", s)
  s <- sub("(\\d)T(\\d)", "\\1 \\2", s)
  s <- sub(" UTC$", "", s)
  s <- sub(" (\\d{1,2})(\\d{2})$", " \\1:\\2", s)      # "2025-01-01 1500" -> 15:00
  out <- as.POSIXct(rep(NA_real_, length(s)), origin = "1970-01-01", tz = "UTC")
  ano_primeiro <- grepl("^\\d{4}[-/]", s)
  dia_primeiro <- grepl("^\\d{1,2}/\\d{1,2}/\\d{4}", s)
  tenta <- function(idx, formatos) {
    for (f in formatos) {
      i <- idx & is.na(out)
      if (!any(i)) break
      out[i] <<- as.POSIXct(s[i], format = f, tz = "UTC")
    }
  }
  tenta(ano_primeiro, c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y/%m/%d %H:%M:%S",
                        "%Y/%m/%d %H:%M", "%Y-%m-%d %H", "%Y-%m-%d"))
  tenta(dia_primeiro, c("%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y"))
  if (anyNA(out[!is.na(s) & nzchar(s)])) {
    ex <- s[is.na(out) & !is.na(s)][1]
    warning("Algumas datas não foram reconhecidas (exemplo: '", ex, "').", call. = FALSE)
  }
  out
}

msg("Saída em: ", DIR_SAIDA)
dir.create(DIR_SAIDA, recursive = TRUE, showWarnings = FALSE)


## ============================================================================
## Etapa A. Observações: leitura, referência temporal e controle de qualidade
## ============================================================================
msg("Etapa A: observações")
dA <- pasta("01_observacoes")

# ---- A.1 Leitores --------------------------------------------------------------------
ler_inmet_oficial <- function(pasta_inmet) {
  arqs <- list.files(pasta_inmet, pattern = "\\.csv$", ignore.case = TRUE,
                     recursive = TRUE, full.names = TRUE)
  arqs <- arqs[grepl("_RN_", toupper(basename(arqs)))]
  if (!length(arqs)) stop("Nenhum CSV do RN (nome com '_RN_') em ", pasta_inmet, call. = FALSE)
  purrr::map_dfr(arqs, function(f) {
    cab <- iconv(readLines(f, n = 8, warn = FALSE), "latin1", "UTF-8")
    campo <- function(rotulo) {
      l <- cab[grepl(rotulo, chave(cab))][1]
      trimws(gsub(";+\\s*$", "", sub("^[^;]*;", "", l)))
    }
    nome <- campo("^ESTA"); codigo <- campo("^CODIGO")
    lat <- num(campo("^LATITUDE")); lon <- num(campo("^LONGITUDE"))
    d <- utils::read.csv2(f, skip = 8, header = TRUE, fileEncoding = "latin1",
                          check.names = FALSE, stringsAsFactors = FALSE,
                          na.strings = c("", "-9999", "-9999,0", "-9999.0"))
    nm <- chave(names(d))
    col <- function(padrao) { j <- grep(padrao, nm); if (length(j)) j[1] else NA_integer_ }
    valor <- function(j) if (is.na(j)) rep(NA_real_, nrow(d)) else num(d[[j]])
    jd <- col("^DATA"); jh <- col("^HORA")
    hh <- sprintf("%04d", suppressWarnings(as.integer(substr(gsub("\\D", "", d[[jh]]), 1, 4))))
    dd <- gsub("-", "/", as.character(d[[jd]]))
    dt <- ifelse(grepl("^\\d{4}/", dd),
                 format(as.POSIXct(paste(dd, hh), format = "%Y/%m/%d %H%M", tz = "UTC"), "%Y-%m-%d %H:%M"),
                 format(as.POSIXct(paste(dd, hh), format = "%d/%m/%Y %H%M", tz = "UTC"), "%Y-%m-%d %H:%M"))
    tibble(estacao = paste0(toupper(nome), " (", codigo, ")"), rede = "INMET",
           lon = lon, lat = lat,
           datahora_utc = as.POSIXct(dt, format = "%Y-%m-%d %H:%M", tz = "UTC"),
           T = valor(col("BULBO SECO")),
           UR = valor(col("UMIDADE RELATIVA DO AR, HORARIA")),
           V = valor(col("VENTO, VELOCIDADE HORARIA")),
           RAD = valor(col("RADIACAO GLOBAL")),
           altura_vento = 10)
  })
}

ler_emparn_bruto <- function(arq, cols) {
  l1 <- readLines(arq, n = 1, warn = FALSE)
  n_pv <- lengths(regmatches(l1, gregexpr(";", l1)))
  n_vg <- lengths(regmatches(l1, gregexpr(",", l1)))
  d <- if (n_pv > n_vg) {
    readr::read_delim(arq, delim = ";", locale = readr::locale(decimal_mark = ","),
                      show_col_types = FALSE, guess_max = 100000)
  } else {
    readr::read_csv(arq, show_col_types = FALSE, guess_max = 100000)
  }
  precisa <- unlist(cols[c("estacao", "lon", "lat", "datahora", "T", "UR", "V")])
  faltam <- setdiff(precisa, names(d))
  if (length(faltam)) {
    stop("Colunas não encontradas no arquivo da EMPARN: ", paste(faltam, collapse = ", "),
         "\nColunas disponíveis: ", paste(names(d), collapse = ", "),
         "\nAjuste COLS_EMPARN no bloco de configuração.", call. = FALSE)
  }
  carimbo <- ler_datahora(do.call(paste, unname(as.list(d[cols$datahora]))))
  rad <- if (!is.na(cols$RAD) && cols$RAD %in% names(d)) num(d[[cols$RAD]]) else NA_real_
  tibble(estacao = as.character(d[[cols$estacao]]), rede = "EMPARN",
         lon = num(d[[cols$lon]]), lat = num(d[[cols$lat]]), carimbo = carimbo,
         T = num(d[[cols$T]]), UR = num(d[[cols$UR]]), V = num(d[[cols$V]]),
         RAD = rad, altura_vento = 2)
}

# ---- A.2 Leitura conforme o modo ----------------------------------------------------
if (MODO_OBS == "brutos") {
  if (BAIXAR_INMET && !dir.exists(DIR_INMET_BRUTO)) {
    dir.create(DIR_INMET_BRUTO, recursive = TRUE, showWarnings = FALSE)
    zipf <- file.path(DIR_INMET_BRUTO, paste0(P$ano, ".zip"))
    utils::download.file(paste0("https://portal.inmet.gov.br/uploads/dadoshistoricos/", P$ano, ".zip"),
                         zipf, mode = "wb")
    utils::unzip(zipf, exdir = DIR_INMET_BRUTO)
  }
  inm <- ler_inmet_oficial(DIR_INMET_BRUTO)
  emp <- ler_emparn_bruto(ARQ_EMPARN_BRUTO, COLS_EMPARN)

  # Verificação da referência temporal (Seção 2.2): hora do máximo médio da radiação
  # e da temperatura em cada leitura do carimbo.
  leituras <- dplyr::bind_rows(
    emp %>% mutate(leitura = "EMPARN, carimbo lido como UTC",
                   hl = as.integer(format(carimbo + P$fuso * 3600, "%H"))),
    emp %>% mutate(leitura = "EMPARN, carimbo lido como hora local",
                   hl = as.integer(format(carimbo, "%H"))),
    inm %>% mutate(leitura = "INMET, UTC convertido para hora local",
                   hl = as.integer(format(datahora_utc + P$fuso * 3600, "%H"))))
  ciclo_ref <- leituras %>%
    mutate(RAD = ifelse(is.finite(RAD) & RAD >= 0, RAD, NA)) %>%
    group_by(leitura, hora_local = hl) %>%
    summarise(radiacao_media = mean(RAD, na.rm = TRUE),
              T_media = mean(ifelse(T >= P$qc$T[1] & T <= P$qc$T[2], T, NA), na.rm = TRUE),
              .groups = "drop")
  checagem <- ciclo_ref %>% group_by(leitura) %>%
    summarise(hora_max_radiacao = if (all(is.na(radiacao_media))) NA_integer_ else
                hora_local[which.max(radiacao_media)],
              hora_max_temperatura = hora_local[which.max(T_media)], .groups = "drop")
  gravar(ciclo_ref, dA, "referencia_temporal_ciclo_medio.csv")
  gravar(checagem, dA, "referencia_temporal_checagem.csv")
  print(as.data.frame(checagem))
  # Esperado (artigo): radiação máxima às 9h se o carimbo da EMPARN for lido como UTC,
  # às 12h se lido como hora local; temperatura máxima às 14h (EMPARN) e 13h (INMET).

  emp <- emp %>% mutate(datahora_utc = if (CARIMBO_EMPARN_EM_HORA_LOCAL)
    carimbo - P$fuso * 3600 else carimbo)
  obs_bruto <- dplyr::bind_rows(
    emp %>% select(estacao, rede, lon, lat, datahora_utc, T, UR, V, altura_vento),
    inm %>% select(estacao, rede, lon, lat, datahora_utc, T, UR, V, altura_vento))
  gravar(obs_bruto %>% mutate(datahora = format(datahora_utc, "%Y-%m-%dT%H:%M:%SZ")) %>%
           select(estacao, rede, lon, lat, datahora, T, UR, V, altura_vento),
         dA, "obs_horarias_padronizadas_2025.csv")
} else if (MODO_OBS == "padronizado") {
  if (!existe(ARQ_OBS_PADRONIZADO)) stop("Arquivo não encontrado: ", ARQ_OBS_PADRONIZADO, call. = FALSE)
  d <- readr::read_csv(ARQ_OBS_PADRONIZADO, show_col_types = FALSE, guess_max = 100000)
  faltam <- setdiff(unname(COLS_OBS), names(d))
  if (length(faltam)) {
    stop("Colunas não encontradas em ", basename(ARQ_OBS_PADRONIZADO), ": ",
         paste(faltam, collapse = ", "), "\nColunas disponíveis: ",
         paste(names(d), collapse = ", "), "\nAjuste COLS_OBS.", call. = FALSE)
  }
  obs_bruto <- tibble(estacao = as.character(d[[COLS_OBS[["estacao"]]]]),
                      rede = toupper(as.character(d[[COLS_OBS[["rede"]]]])),
                      lon = num(d[[COLS_OBS[["lon"]]]]), lat = num(d[[COLS_OBS[["lat"]]]]),
                      datahora_utc = ler_datahora(d[[COLS_OBS[["datahora"]]]]),
                      T = num(d[[COLS_OBS[["T"]]]]), UR = num(d[[COLS_OBS[["UR"]]]]),
                      V = num(d[[COLS_OBS[["V"]]]]))
  obs_bruto$altura_vento <- if ("altura_vento" %in% names(d)) num(d$altura_vento) else
    ifelse(obs_bruto$rede == "INMET", 10, 2)
  if (V_JA_A_2M) {   # devolve o vento do INMET à altura nativa (10 m)
    obs_bruto <- obs_bruto %>%
      mutate(V = ifelse(altura_vento != 2, V / fator_altura(altura_vento), V))
  }
} else stop("MODO_OBS deve ser 'padronizado' ou 'brutos'.", call. = FALSE)

if (any(is.na(obs_bruto$datahora_utc))) {
  warning(sum(is.na(obs_bruto$datahora_utc)), " linhas sem data e hora reconhecidas foram descartadas.",
          call. = FALSE)
  obs_bruto <- obs_bruto %>% filter(!is.na(datahora_utc))
}

# ---- A.3 Catálogo das estações ------------------------------------------------------
catalogo <- obs_bruto %>% group_by(estacao, rede) %>%
  summarise(lon = dplyr::first(stats::na.omit(lon)), lat = dplyr::first(stats::na.omit(lat)),
            .groups = "drop")
if (existe(ARQ_CATALOGO)) {
  cat_ext <- readr::read_csv(ARQ_CATALOGO, show_col_types = FALSE) %>%
    transmute(estacao = as.character(estacao), rede = toupper(rede), lon = num(lon), lat = num(lat))
  catalogo <- dplyr::bind_rows(catalogo, cat_ext %>% filter(!chave(estacao) %in% chave(catalogo$estacao)))
}
catalogo <- ordenar_estacoes(catalogo)
gravar(catalogo, dA, "catalogo_estacoes.csv")
msg("   ", nrow(catalogo), " estações no catálogo (", sum(catalogo$rede == "EMPARN"), " EMPARN, ",
    sum(catalogo$rede == "INMET"), " INMET)")

# ---- A.4 Hora local, ano e controle de qualidade --------------------------------------
dentro <- function(x, lim) ifelse(!is.na(x) & x >= lim[1] & x <= lim[2], x, NA_real_)
obs <- obs_bruto %>%
  distinct(estacao, datahora_utc, .keep_all = TRUE) %>%
  mutate(datahora_local = datahora_utc + P$fuso * 3600,
         data = as.Date(datahora_local), hora = as.integer(format(datahora_local, "%H")),
         ano = as.integer(format(datahora_local, "%Y"))) %>%
  filter(ano == P$ano) %>%
  mutate(T = dentro(T, P$qc$T), UR = dentro(UR, P$qc$UR), V = dentro(V, P$qc$V),
         valido = !is.na(T) & !is.na(UR) & !is.na(V))

resumo_qc <- obs %>% group_by(estacao, rede) %>%
  summarise(registros = n(), T_validos = sum(!is.na(T)), UR_validos = sum(!is.na(UR)),
            V_validos = sum(!is.na(V)), validos_conjuntos = sum(valido), .groups = "drop") %>%
  right_join(catalogo %>% select(estacao, rede), by = c("estacao", "rede")) %>%
  mutate(across(registros:validos_conjuntos, ~ dplyr::coalesce(.x, 0L))) %>%
  ordenar_estacoes()
gravar(resumo_qc, dA, "controle_qualidade_por_estacao.csv")


## ============================================================================
## Etapa B. Índices por observação e médias por estação e horário
## ============================================================================
msg("Etapa B: índices e médias por estação")
dB <- pasta("02_estacoes")

# Os índices são calculados em cada observação e só depois promediados (Seção 2.2).
val <- obs %>% filter(valido) %>%
  mutate(V2 = V * fator_altura(altura_vento),
         IC = ic_c(T, UR), NET = net_c(T, UR, V2), NET0 = net_c(T, UR, 0),
         dif = IC - NET, p_tur = IC - NET0, p_vento = NET0 - NET)

medias_hora <- function(h) {
  m <- val %>% filter(hora == h) %>% group_by(estacao, rede) %>%
    summarise(n = n(), T = mean(T), UR = mean(UR), V = mean(V2), IC = mean(IC),
              NET = mean(NET), NET0 = mean(NET0), .groups = "drop")
  catalogo %>% left_join(m, by = c("estacao", "rede")) %>%
    mutate(n = dplyr::coalesce(n, 0L), horario = h, retida = n >= P$min_registros,
           situacao = dplyr::case_when(retida ~ "Retida na análise",
                                       n > 0 ~ paste0("Abaixo de ", P$min_registros, " registros"),
                                       TRUE ~ "Sem registro conjunto das três variáveis")) %>%
    ordenar_estacoes()
}
est_h <- purrr::set_names(purrr::map(P$horas, medias_hora), as.character(P$horas))
for (h in names(est_h)) gravar(est_h[[h]], dB, sprintf("estacoes_h%02d.csv", as.integer(h)))
retidas <- purrr::map(est_h, ~ .x$estacao[.x$retida])

resumo_rede <- purrr::map_dfr(est_h, function(e) tibble(
  horario = e$horario[1], estacoes_retidas = sum(e$retida),
  retidas_EMPARN = sum(e$retida & e$rede == "EMPARN"), retidas_INMET = sum(e$retida & e$rede == "INMET"),
  registros_retidas = sum(e$n[e$retida]), abaixo_do_limiar = sum(!e$retida & e$n > 0),
  sem_registro = sum(e$n == 0)))
gravar(resumo_rede, dB, "resumo_rede_por_horario.csv")
print(as.data.frame(resumo_rede))

# Comparação opcional com as médias originais do artigo
comparar_original <- function(novo, arq, chaves = "estacao", colunas) {
  if (!existe(arq)) return(invisible(NULL))
  orig <- readr::read_csv(arq, show_col_types = FALSE)
  colunas <- intersect(colunas, intersect(names(orig), names(novo)))
  j <- dplyr::inner_join(novo %>% mutate(k = chave(estacao)) %>% select(k, all_of(colunas)),
                         orig %>% mutate(k = chave(estacao)) %>% select(k, all_of(colunas)),
                         by = "k", suffix = c("_novo", "_orig"))
  purrr::map_dfr(colunas, function(cc) tibble(
    variavel = cc, estacoes = nrow(j),
    maior_diferenca_absoluta = max(abs(num(j[[paste0(cc, "_novo")]]) - num(j[[paste0(cc, "_orig")]])),
                                   na.rm = TRUE)))
}
cmp15 <- comparar_original(est_h[["15"]], file.path(DIR_ORIGINAL, "estacoes_h15.csv"),
                           colunas = c("n", "T", "UR", "V", "IC", "NET", "NET0"))
if (!is.null(cmp15)) {
  gravar(cmp15, dB, "comparacao_com_original_h15.csv")
  msg("   comparação com estacoes_h15.csv original:")
  print(as.data.frame(cmp15))
}


## ============================================================================
## Etapa C. Limites, grade de 500 m e covariáveis
## ============================================================================
msg("Etapa C: limites, grade e covariáveis")
dC <- pasta("03_covariaveis")

ler_vetor <- function(arq) if (existe(arq)) sf::st_read(arq, quiet = TRUE) else NULL
remover_buracos <- function(g) {
  partes <- sf::st_cast(sf::st_geometry(g), "POLYGON")
  sf::st_union(sf::st_sfc(lapply(partes, function(p) sf::st_polygon(list(p[[1]]))),
                          crs = sf::st_crs(g)))
}

# ---- C.1 Grade -------------------------------------------------------------------------
rn_lim <- ler_vetor(ARQ_LIMITE_RN)
mun_lim <- ler_vetor(ARQ_MUNICIPIOS)
costa_lim <- ler_vetor(ARQ_LINHA_COSTA)
if ((is.null(rn_lim) || is.null(mun_lim) || is.null(costa_lim)) && !tem("geobr")) {
  stop("Instale o pacote geobr ou informe ARQ_LIMITE_RN, ARQ_MUNICIPIOS e ARQ_LINHA_COSTA.",
       call. = FALSE)
}
if (is.null(rn_lim)) rn_lim <- geobr::read_state(code_state = "RN", year = ANO_MALHA, showProgress = FALSE)
if (is.null(mun_lim)) mun_lim <- geobr::read_municipality(code_muni = "RN", year = ANO_MALHA,
                                                         showProgress = FALSE)

if (existe(ARQ_GRADE_MODELO)) {
  grade <- terra::rast(terra::rast(ARQ_GRADE_MODELO))      # só a geometria
  msg("   grade copiada de ", basename(ARQ_GRADE_MODELO))
} else {
  rn_tmp <- sf::st_transform(sf::st_make_valid(rn_lim), sf::st_crs(P$crs))
  grade <- terra::rast(terra::ext(terra::vect(rn_tmp)), nrows = P$nrow, ncols = P$ncol, crs = P$crs)
}
CRS_GRADE <- sf::st_crs(terra::crs(grade))

rn <- sf::st_transform(sf::st_make_valid(rn_lim), CRS_GRADE) %>% sf::st_union()
rn_cheio <- remover_buracos(rn)                 # limite sem as lagoas e salinas
cod_col <- intersect(c("code_muni", "CD_MUN", "cod_ibge", "codigo", "CD_GEOCMU"), names(mun_lim))[1]
nom_col <- intersect(c("name_muni", "NM_MUN", "nome", "municipio", "NM_MUNICIP"), names(mun_lim))[1]
if (is.na(cod_col) || is.na(nom_col)) stop("Não encontrei código e nome nos municípios.", call. = FALSE)
mun <- sf::st_transform(sf::st_make_valid(mun_lim), CRS_GRADE)
mun <- sf::st_sf(code_muni = as.numeric(as.character(mun[[cod_col]])),
                 name_muni = as.character(mun[[nom_col]]), geometry = sf::st_geometry(mun))

if (is.null(costa_lim)) {
  viz <- dplyr::bind_rows(
    geobr::read_state(code_state = "CE", year = ANO_MALHA, showProgress = FALSE),
    geobr::read_state(code_state = "PB", year = ANO_MALHA, showProgress = FALSE)) %>%
    sf::st_transform(CRS_GRADE) %>% sf::st_make_valid() %>% sf::st_union()
  costa <- sf::st_difference(sf::st_boundary(rn_cheio), sf::st_buffer(viz, 1000))
} else {
  costa <- sf::st_union(sf::st_transform(costa_lim, CRS_GRADE))
}
sf::st_write(sf::st_sf(geometry = costa), file.path(dC, "linha_costa.gpkg"), quiet = TRUE,
             delete_dsn = TRUE)

mascara <- terra::rasterize(terra::vect(sf::st_sf(geometry = rn)), grade, field = 1)
mascara_exib <- terra::rasterize(terra::vect(sf::st_sf(geometry = rn_cheio)), grade, field = 1)
names(mascara) <- "mascara"; names(mascara_exib) <- "mascara_exib"
info_grade <- tibble(
  colunas = terra::ncol(grade), linhas = terra::nrow(grade),
  resolucao_x_m = terra::res(grade)[1], resolucao_y_m = terra::res(grade)[2],
  xmin = terra::xmin(grade), xmax = terra::xmax(grade), ymin = terra::ymin(grade), ymax = terra::ymax(grade),
  celulas_analise = sum(!is.na(terra::values(mascara))),
  celulas_exibicao = sum(!is.na(terra::values(mascara_exib))))
gravar(info_grade, dC, "grade_500m.csv")
print(as.data.frame(info_grade))

# ---- C.2 Covariáveis fixas ---------------------------------------------------------------
if (BAIXAR_LCZ_COM_LCZ4R && !existe(ARQ_LCZ)) {
  if (!tem("LCZ4r")) stop("Instale o pacote LCZ4r ou forneça lcz_rn.tif.", call. = FALSE)
  lcz_bruto <- LCZ4r::lcz_get_map(roi = sf::st_sf(geometry = sf::st_transform(rn_cheio, 4326)))
  terra::writeRaster(lcz_bruto, ARQ_LCZ, overwrite = TRUE)
}
x_r <- terra::init(grade, "x"); names(x_r) <- "x"
y_r <- terra::init(grade, "y"); names(y_r) <- "y"
if (!USAR_PILHA_ORIGINAL) {
  for (f in c(ARQ_SRTM, ARQ_NDVI, ARQ_LCZ, ARQ_ERA5)) {
    if (!existe(f)) stop("Covariável não encontrada: ", f, call. = FALSE)
  }
  # Altitude e NDVI: interpolação bilinear para a grade; LCZ: classe mais frequente.
  alt_r <- terra::project(terra::rast(ARQ_SRTM)[[1]], grade, method = P$reamostragem_continua)
  ndvi_r <- terra::project(terra::rast(ARQ_NDVI)[[1]], grade, method = P$reamostragem_continua)
  lcz_r <- terra::project(terra::rast(ARQ_LCZ)[[1]], grade, method = P$reamostragem_lcz)
  names(alt_r) <- "altitude"; names(ndvi_r) <- "ndvi"; names(lcz_r) <- "lcz"
  # Distância euclidiana à linha de costa (km), medida do centro de cada célula.
  pts_grade <- terra::vect(terra::xyFromCell(grade, seq_len(terra::ncell(grade))), crs = terra::crs(grade))
  d_costa <- terra::distance(pts_grade, terra::aggregate(terra::vect(sf::st_sf(geometry = costa))))
  dist_r <- terra::rast(grade)
  terra::values(dist_r) <- as.numeric(d_costa) / 1000
  names(dist_r) <- "dist_costa"
  rm(pts_grade, d_costa)
  fixas <- c(x_r, y_r, alt_r, dist_r, ndvi_r, lcz_r)
}

# ---- C.3 ERA5-Land por horário e pilhas -------------------------------------------------
pilha_hora <- function(h) {
  if (USAR_PILHA_ORIGINAL) {
    f <- sprintf(ARQ_PILHA_ORIGINAL, (h - P$fuso) %% 24)          # 00h -> 03 UTC etc.
    if (!existe(f)) stop("Pilha original não encontrada: ", f, call. = FALSE)
    r <- terra::rast(f)
    precisa <- c("altitude", "dist_costa", "ndvi", "lcz", unname(COL_ERA5))
    falt <- setdiff(precisa, names(r))
    if (length(falt)) stop("Bandas ausentes em ", basename(f), ": ", paste(falt, collapse = ", "),
                           "\nBandas disponíveis: ", paste(names(r), collapse = ", "), call. = FALSE)
    r <- r[[precisa]]
    if (!terra::compareGeom(r, grade, stopOnError = FALSE)) r <- terra::project(r, grade, method = "near")
    return(terra::mask(c(x_r, y_r, r), mascara))
  }
  r <- terra::rast(ARQ_ERA5[[as.character(h)]])
  falt <- setdiff(unname(COL_ERA5), names(r))
  if (length(falt)) stop("Bandas ausentes em ", basename(ARQ_ERA5[[as.character(h)]]), ": ",
                         paste(falt, collapse = ", "), call. = FALSE)
  e <- terra::project(r[[unname(COL_ERA5)]], grade, method = "bilinear")
  terra::mask(c(fixas, e), mascara)
}
pilhas <- purrr::set_names(purrr::map(P$horas, pilha_hora), as.character(P$horas))
for (h in names(pilhas)) {
  terra::writeRaster(pilhas[[h]], file.path(dC, sprintf("covariaveis_500m_h%02d.tif", as.integer(h))),
                     overwrite = TRUE)
}
terra::writeRaster(c(mascara, mascara_exib), file.path(dC, "mascaras_500m.tif"), overwrite = TRUE)

# ---- C.4 Covariáveis nas estações ---------------------------------------------------------
# Valor da célula que contém a estação; se faltar, média das células válidas num raio
# de 1, 3, 6 e 12 km (classe mais frequente para a LCZ). x e y são as coordenadas exatas.
pontos_estacoes <- function(est) {
  sf::st_as_sf(est, coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    sf::st_transform(CRS_GRADE)
}
cov_estacoes <- function(pilha, est) {
  pts <- pontos_estacoes(est)
  xy <- sf::st_coordinates(pts)
  v <- terra::extract(pilha, terra::vect(pts), ID = FALSE)
  for (cc in setdiff(names(pilha), c("x", "y"))) {
    for (i in which(is.na(v[[cc]]))) {
      for (r in P$raios_busca) {
        vals <- terra::extract(pilha[[cc]], terra::buffer(terra::vect(pts[i, ]), r), ID = FALSE)[[1]]
        vals <- vals[!is.na(vals)]
        if (length(vals)) { v[[cc]][i] <- if (cc == "lcz") moda(vals) else mean(vals); break }
      }
    }
  }
  v$x <- xy[, 1]; v$y <- xy[, 2]
  dplyr::bind_cols(est, tibble::as_tibble(v))
}
est_cov <- purrr::imap(est_h, function(e, h) cov_estacoes(pilhas[[h]], e %>% filter(retida)))
for (h in names(est_cov)) gravar(est_cov[[h]], dC, sprintf("estacoes_cov_h%02d.csv", as.integer(h)))
cmpc <- comparar_original(est_cov[["15"]], file.path(DIR_ORIGINAL, "estacoes_cov_h15.csv"),
                          colunas = c("altitude", "dist_costa", "ndvi", "lcz", unname(COL_ERA5), "x", "y"))
if (!is.null(cmpc)) {
  gravar(cmpc, dC, "comparacao_covariaveis_com_original_h15.csv")
  msg("   comparação com estacoes_cov_h15.csv original:")
  print(as.data.frame(cmpc))
}


## ============================================================================
## Etapa D. Validação cruzada (M1 a M15, R1 e R2) nos três horários
## ============================================================================
msg("Etapa D: validação cruzada leave-one-out")
dD <- pasta("04_validacao")

como_sp <- function(d) {
  d <- as.data.frame(d)
  sp::coordinates(d) <- ~ x + y
  d
}
.krige <- function(form, tr, nd) {
  trs <- como_sp(tr); nds <- como_sp(nd)
  vm <- suppressWarnings(automap::autofitVariogram(form, trs)$var_model)
  suppressWarnings(gstat::krige(form, trs, nds, model = vm, debug.level = 0)$var1.pred)
}
.idw <- function(tr, nd) gstat::idw(valor ~ 1, como_sp(tr), como_sp(nd), idp = 2, debug.level = 0)$var1.pred
.lm <- function(form, tr, nd) unname(stats::predict(stats::lm(form, data = tr), newdata = nd))
.lm_lcz <- function(tr, nd) {
  tr$lcz <- factor(tr$lcz)
  nd$lcz <- factor(nd$lcz, levels = levels(tr$lcz))   # classe ausente no treino -> NA
  fit <- stats::lm(valor ~ x + y + altitude + lcz, data = tr)
  unname(suppressWarnings(stats::predict(fit, newdata = nd)))
}
.gam <- function(form, tr, nd) {
  fit <- mgcv::gam(form, data = tr, method = "REML")
  as.numeric(stats::predict(fit, newdata = nd))
}
.rf <- function(tr, nd) {
  fit <- ranger::ranger(valor ~ x + y + altitude + era5,
                        data = as.data.frame(tr[, c("valor", "x", "y", "altitude", "era5")]),
                        num.trees = 500, min.node.size = 2, seed = 1)
  stats::predict(fit, data = as.data.frame(nd[, c("x", "y", "altitude", "era5")]))$predictions
}

# Ajusta o método 'codigo' com 'tr' (colunas valor e covariáveis) e prediz 'nd'.
prever <- function(codigo, tr, nd) {
  nec <- NECESSARIAS[[codigo]]
  tr <- tr[stats::complete.cases(tr[, c("valor", nec)]), , drop = FALSE]
  ok <- stats::complete.cases(nd[, nec, drop = FALSE])
  out <- rep(NA_real_, nrow(nd))
  if (!any(ok) || nrow(tr) < 3) return(out)
  ndd <- nd[ok, , drop = FALSE]
  p <- tryCatch(switch(codigo,
    IDW                = .idw(tr, ndd),
    OK                 = .krige(valor ~ 1, tr, ndd),
    KDE_alt            = .krige(valor ~ altitude, tr, ndd),
    KDE_alt_ndvi       = .krige(valor ~ altitude + ndvi, tr, ndd),
    LM_coords          = .lm(valor ~ x + y, tr, ndd),
    LM_coords_alt      = .lm(valor ~ x + y + altitude, tr, ndd),
    LM_coords_alt_lcz  = .lm_lcz(tr, ndd),
    GAM_coords         = .gam(valor ~ te(x, y, k = 3), tr, ndd),
    GAM_coords_alt     = .gam(valor ~ te(x, y, k = 3) + s(altitude, k = 3), tr, ndd),
    LM_alt_era5        = .lm(valor ~ altitude + era5, tr, ndd),
    LM_coords_alt_era5 = .lm(valor ~ x + y + altitude + era5, tr, ndd),
    KDE_alt_era5       = .krige(valor ~ altitude + era5, tr, ndd),
    RF_era5            = .rf(tr, ndd),
    LM_alt_era5_dist   = .lm(valor ~ altitude + era5 + dist_costa, tr, ndd),
    KDE_alt_era5_dist  = .krige(valor ~ altitude + era5 + dist_costa, tr, ndd),
    stop("Método desconhecido: ", codigo)),
    error = function(e) rep(NA_real_, nrow(ndd)))
  if (length(p) != nrow(ndd)) p <- rep(NA_real_, nrow(ndd))
  out[ok] <- as.numeric(p)
  out
}

dados_variavel <- function(est, v) {
  est %>% transmute(estacao, valor = .data[[v]], x, y, altitude, ndvi, lcz, dist_costa,
                    era5 = .data[[COL_ERA5[[v]]]]) %>%
    filter(!is.na(valor))
}
validacao_cruzada <- function(est, rotulo = "") {
  purrr::map_dfr(VARIAVEIS, function(v) {
    msg("   ", rotulo, " ", v)
    d <- dados_variavel(est, v)
    purrr::map_dfr(seq_len(nrow(d)), function(i) {
      tr <- d[-i, ]; te <- d[i, ]
      tibble(variavel = v, estacao = te$estacao, obs = te$valor, codigo = METODOS$codigo,
             pred = vapply(METODOS$codigo, function(m) prever(m, tr, te), numeric(1)))
    })
  }) %>% left_join(METODOS %>% select(codigo, metodo = id), by = "codigo") %>%
    select(variavel, metodo, estacao, obs, pred)
}
# R1: ERA5-Land no ponto; R2: R1 + viés médio (obs - ERA5) das DEMAIS estações.
referencias_reanalise <- function(est) {
  purrr::map_dfr(VARIAVEIS, function(v) {
    d <- est %>% transmute(estacao, obs = .data[[v]], era5 = .data[[COL_ERA5[[v]]]]) %>%
      filter(!is.na(obs), !is.na(era5))
    vies_demais <- (sum(d$obs - d$era5) - (d$obs - d$era5)) / (nrow(d) - 1)
    dplyr::bind_rows(
      tibble(variavel = v, metodo = "R1", estacao = d$estacao, obs = d$obs, pred = d$era5),
      tibble(variavel = v, metodo = "R2", estacao = d$estacao, obs = d$obs, pred = d$era5 + vies_demais))
  })
}
metricas <- function(cv) {
  cv %>% group_by(variavel, horario, metodo) %>%
    summarise(n = sum(!is.na(pred)), n_falhas = sum(is.na(pred)),
              R2 = { ok <- !is.na(pred); 1 - sum((pred[ok] - obs[ok])^2) / sum((obs[ok] - mean(obs[ok]))^2) },
              RMSE = rmse(pred - obs), MAE = mean(abs(pred - obs), na.rm = TRUE),
              vies = mean(pred - obs, na.rm = TRUE), .groups = "drop")
}

cv <- em_cache(file.path(dD, "validacao_cruzada.rds"), {
  purrr::map_dfr(P$horas, function(h) {
    e <- est_cov[[as.character(h)]]
    dplyr::bind_rows(validacao_cruzada(e, sprintf("%02dh", h)), referencias_reanalise(e)) %>%
      mutate(horario = h)
  }) %>% select(variavel, horario, metodo, estacao, obs, pred)
}, ass = do.call(assinatura, unname(est_cov)))
gravar(cv, dD, "loocv_longo.csv")

met <- metricas(cv) %>% left_join(ROTULOS %>% select(metodo = id, codigo, familia, rotulo), by = "metodo")
gravar(met, dD, "metricas_loocv.csv")

# Vencedor: menor RMSE entre M1 e M15 que predisseram todas as estações.
venc <- met %>% filter(metodo %in% METODOS$id, n_falhas == 0) %>%
  group_by(variavel, horario) %>% slice_min(RMSE, n = 1, with_ties = FALSE) %>% ungroup() %>%
  left_join(met %>% filter(metodo %in% METODOS$id) %>% group_by(variavel, horario) %>%
              summarise(todos_R2_negativos = all(R2 < 0, na.rm = TRUE), .groups = "drop"),
            by = c("variavel", "horario")) %>%
  mutate(sem_poder_preditivo = todos_R2_negativos)
gravar(venc, dD, "vencedores.csv")
print(as.data.frame(venc %>% select(variavel, horario, metodo, codigo, R2, RMSE, sem_poder_preditivo)))
codigo_vencedor <- function(v, h) venc$codigo[venc$variavel == v & venc$horario == h]
metodo_vencedor <- function(v, h) venc$metodo[venc$variavel == v & venc$horario == h]


## ============================================================================
## Etapa E. Incerteza do RMSE: bootstrap e testes pareados
## ============================================================================
msg("Etapa E: bootstrap e testes pareados")
dE <- pasta("05_incerteza")

inc <- em_cache(file.path(dE, "bootstrap.rds"), {
  set.seed(P$semente_bootstrap)
  dados_cv <- cv %>% mutate(erro = pred - obs)
  B <- P$B
  ic_rmse <- dados_cv %>%
    group_by(variavel, horario, metodo) %>%
    group_modify(~ {
      e <- .x$erro; n <- length(e)
      bs <- replicate(B, rmse(sample(e, n, replace = TRUE)))
      tibble(n = n, rmse = rmse(e), ic95_inf = stats::quantile(bs, .025, na.rm = TRUE),
             ic95_sup = stats::quantile(bs, .975, na.rm = TRUE))
    }) %>% ungroup()
  # Cada método contra R2, às 15h: bootstrap pareado da diferença de RMSE e Wilcoxon
  # sobre os erros absolutos.
  wide <- dados_cv %>% filter(horario == P$hora_mapa) %>%
    select(variavel, metodo, estacao, erro) %>%
    pivot_wider(names_from = metodo, values_from = erro)
  testes <- purrr::map_dfr(unique(wide$variavel), function(v) {
    w <- wide %>% filter(variavel == v)
    metodos <- setdiff(names(w), c("variavel", "estacao", "R2"))
    purrr::map_dfr(metodos, function(m) {
      ok <- stats::complete.cases(w[[m]], w[["R2"]])
      a <- w[[m]][ok]; b <- w[["R2"]][ok]; n <- length(a)
      if (n < 5) return(NULL)
      dif <- replicate(B, { i <- sample.int(n, n, TRUE); rmse(a[i]) - rmse(b[i]) })
      tibble(variavel = v, metodo = m, n = n, dif_rmse = rmse(a) - rmse(b),
             ic95_inf = stats::quantile(dif, .025), ic95_sup = stats::quantile(dif, .975),
             p_wilcoxon = tryCatch(stats::wilcox.test(abs(a), abs(b), paired = TRUE,
                                                      exact = FALSE)$p.value,
                                   error = function(e) NA_real_))
    })
  })
  list(ic_rmse = ic_rmse, testes = testes)
}, ass = assinatura(cv))
gravar(inc$ic_rmse, dE, "tabelaS2_rmse_ic95.csv")
gravar(inc$testes, dE, "testes_pareados_vs_R2.csv")


## ============================================================================
## Etapa F. Diagnóstico da krigagem (IC às 15h; Tabela S1)
## ============================================================================
msg("Etapa F: diagnóstico da krigagem")
dF <- pasta("06_krigagem")

kg <- em_cache(file.path(dF, "krigagem.rds"), {
  est <- as.data.frame(est_cov[[as.character(P$hora_mapa)]] %>% transmute(estacao, x, y, valor = IC))
  sp_all <- est; sp::coordinates(sp_all) <- ~ x + y
  extrai <- function(v) {
    if (is.null(v)) return(tibble(modelo = NA_character_, pepita = NA_real_,
                                  patamar_parcial = NA_real_, alcance_m = NA_real_, kappa = NA_real_))
    tibble(modelo = as.character(v$model[2]), pepita = v$psill[1], patamar_parcial = v$psill[2],
           alcance_m = v$range[2], kappa = v$kappa[2])
  }
  # (a) ajuste automático a cada rodada (o que os métodos M2 a M4, M12 e M15 fazem);
  # (b) variograma esférico com pepita ajustado uma vez com todas as estações e fixo;
  # (d) máxima verossimilhança restrita com todas as estações.
  v_all <- gstat::variogram(valor ~ 1, sp_all, cutoff = 200000, width = 200000 / 8)
  aviso_fixo <- NA_character_      # o artigo relata que este ajuste não converge
  v_fix <- withCallingHandlers(
    tryCatch(gstat::fit.variogram(v_all, gstat::vgm(psill = stats::var(est$valor), "Sph",
                                                    range = 80000, nugget = 0.1 * stats::var(est$valor))),
             error = function(e) NULL),
    warning = function(w) { aviso_fixo <<- conditionMessage(w); invokeRestart("muffleWarning") })
  v_reml <- tryCatch(suppressWarnings(gstat::fit.variogram.reml(
    valor ~ 1, sp_all, model = gstat::vgm(stats::var(est$valor), "Sph", 80000, 0.1 * stats::var(est$valor)))),
    error = function(e) NULL)
  res <- purrr::map_dfr(seq_len(nrow(est)), function(i) {
    tr <- sp_all[-i, ]; te <- sp_all[i, ]
    va <- tryCatch(suppressWarnings(automap::autofitVariogram(valor ~ 1, tr)$var_model),
                   error = function(e) NULL)
    pred <- function(mod) if (is.null(mod)) NA_real_ else
      tryCatch(suppressWarnings(gstat::krige(valor ~ 1, tr, te, model = mod, debug.level = 0)$var1.pred),
               error = function(e) NA_real_)
    dplyr::bind_rows(
      extrai(va)     %>% mutate(variante = "a_auto_por_rodada", pred = pred(va)),
      extrai(v_fix)  %>% mutate(variante = "b_fixo_todas", pred = pred(v_fix)),
      extrai(v_reml) %>% mutate(variante = "d_REML_todas", pred = pred(v_reml))) %>%
      mutate(estacao = est$estacao[i], obs = est$valor[i])
  })
  list(parametros = res, n_pares = nrow(est) * (nrow(est) - 1) / 2, aviso_fixo = aviso_fixo)
}, ass = assinatura(est_cov[[as.character(P$hora_mapa)]]))
gravar(kg$parametros, dF, "krigagem_parametros_por_rodada.csv")
faixa_obs <- range(kg$parametros$obs)
kg_resumo <- kg$parametros %>% group_by(variante) %>%
  summarise(rmse = rmse(pred - obs), n_predicoes = sum(!is.na(pred)), n_falhas = sum(is.na(pred)),
            modelos = paste(names(table(modelo)), table(modelo), sep = ": ", collapse = "; "),
            pepita_nula = sum(pepita == 0, na.rm = TRUE),
            alcance_mediano_km = stats::median(alcance_m, na.rm = TRUE) / 1000,
            alcance_min_km = min(alcance_m, na.rm = TRUE) / 1000,
            alcance_max_km = max(alcance_m, na.rm = TRUE) / 1000,
            fora_da_faixa_mais_1C = sum(pred < faixa_obs[1] - 1 | pred > faixa_obs[2] + 1, na.rm = TRUE),
            pred_min = suppressWarnings(min(pred, na.rm = TRUE)),
            pred_max = suppressWarnings(max(pred, na.rm = TRUE)), .groups = "drop") %>%
  mutate(faixa_observada = paste(round(faixa_obs, 1), collapse = " a "), pares_de_estacoes = kg$n_pares,
         aviso_do_ajuste = ifelse(variante == "b_fixo_todas", kg$aviso_fixo, NA_character_))
gravar(kg_resumo, dF, "krigagem_resumo_por_variante.csv")
print(as.data.frame(kg_resumo %>% select(variante, rmse, n_falhas, modelos, pepita_nula, alcance_mediano_km)))


## ============================================================================
## Etapa G. Decomposição IC - NET, classes de vento, ciclo diário e sensibilidade
## ============================================================================
msg("Etapa G: decomposição da diferença entre IC e NET")
dG <- pasta("07_decomposicao")
ret15 <- retidas[[as.character(P$hora_mapa)]]
v15 <- val %>% filter(hora == P$hora_mapa, estacao %in% ret15)
v_ret <- val %>% filter(estacao %in% ret15)               # todas as horas, 16 estações

# ---- G.1 Por estação às 15h (Figura 3b) ------------------------------------------------
dec_est <- v15 %>% group_by(estacao, rede) %>%
  summarise(n = n(), across(c(T, UR, V2, IC, NET, NET0, dif, p_tur, p_vento), mean), .groups = "drop") %>%
  left_join(est_cov[[as.character(P$hora_mapa)]] %>% select(estacao, dist_costa), by = "estacao") %>%
  mutate(pct_vento = 100 * p_vento / dif) %>% ordenar_estacoes()
gravar(dec_est, dG, "decomposicao_por_estacao_15h.csv")

correlacoes <- function(e, o) tibble(
  estacoes = nrow(e), observacoes = nrow(o),
  dif_media = mean(e$dif), parcela_vento = mean(e$p_vento), parcela_TUR = mean(e$p_tur),
  pct_vento = 100 * mean(e$p_vento) / mean(e$dif),
  r_V_dif_estacoes = stats::cor(e$V2, e$dif), r_V_dif_horarias = stats::cor(o$V2, o$dif),
  r_pvento_V = stats::cor(e$p_vento, e$V2), r_pvento_T = stats::cor(e$p_vento, e$T),
  r_ptur_T = stats::cor(e$p_tur, e$T), r_ptur_V = stats::cor(e$p_tur, e$V2))
dec_resumo <- dplyr::bind_rows(
  correlacoes(dec_est, v15) %>% mutate(conjunto = "16 estações retidas, 15h"),
  correlacoes(dec_est %>% filter(!eh_jundiai(estacao)), v15 %>% filter(!eh_jundiai(estacao))) %>%
    mutate(conjunto = "sem Jundiaí, 15h")) %>% relocate(conjunto)
gravar(dec_resumo, dG, "decomposicao_resumo_15h.csv")

# ---- G.2 Todas as horas e calmaria ---------------------------------------------------------
todas <- dplyr::bind_rows(
  val %>% summarise(observacoes = n(), r_V_dif = stats::cor(V2, dif)) %>% mutate(conjunto = "todas as estações"),
  v_ret %>% summarise(observacoes = n(), r_V_dif = stats::cor(V2, dif)) %>% mutate(conjunto = "16 estações retidas"))
calma <- dplyr::bind_rows(
  val %>% filter(V2 < 0.5) %>% mutate(conjunto = "todas as estações, todas as horas"),
  v_ret %>% filter(V2 < 0.5) %>% mutate(conjunto = "16 estações retidas, todas as horas"),
  v15 %>% filter(V2 < 0.5) %>% mutate(conjunto = "16 estações retidas, 15h"),
  v15 %>% filter(V2 < 0.5, !eh_jundiai(estacao)) %>% mutate(conjunto = "15h, sem Jundiaí")) %>%
  group_by(conjunto) %>%
  summarise(observacoes = n(), dif_media = mean(dif),
            pct_noturnas = 100 * mean(hora >= 18 | hora < 6),
            pct_jundiai = 100 * mean(eh_jundiai(estacao)), .groups = "drop")
gravar(todas, dG, "correlacao_todas_as_horas.csv")
gravar(calma, dG, "calmaria.csv")

# ---- G.3 Classes de vento às 15h (Tabela 4) ----------------------------------------------
rot_classes <- c("abaixo de 0,5", "0,5 a 1,5", "1,5 a 3", "3 a 4", "4 a 5", "5 ou mais")
classes <- v15 %>%
  mutate(classe = cut(V2, P$classes_vento, labels = rot_classes, right = FALSE)) %>%
  group_by(classe) %>%
  summarise(n = n(), V_media = mean(V2), T = mean(T), UR = mean(UR), IC = mean(IC), NET = mean(NET),
            dif = mean(dif), p_tur = mean(p_tur), p_vento = mean(p_vento),
            n_jundiai = sum(eh_jundiai(estacao)), .groups = "drop") %>%
  mutate(pct = 100 * n / sum(n))
gravar(classes, dG, "classes_vento_15h.csv")

# ---- G.4 Ciclo diário (Figura 4): média das estações em cada hora ----------------------------
ciclo <- v_ret %>% group_by(estacao, hora) %>%
  summarise(across(c(T, V2, dif, p_tur, p_vento), mean), .groups = "drop") %>%
  group_by(hora) %>%
  summarise(estacoes = n(), across(c(T, V2, dif, p_tur, p_vento), mean), .groups = "drop") %>%
  mutate(pct_vento = 100 * p_vento / dif)
gravar(ciclo, dG, "ciclo_diario.csv")

# ---- G.5 Sensibilidade (Seção 2.7) ----------------------------------------------------------
decompor <- function(df, v_ref = 0, z0 = P$z0, excluir_jundiai = FALSE) {
  if (excluir_jundiai) df <- df %>% filter(!eh_jundiai(estacao))
  df <- df %>% mutate(V2s = V * fator_altura(altura_vento, z0), ICs = ic_c(T, UR),
                      NETs = net_c(T, UR, V2s), NET0s = net_c(T, UR, v_ref),
                      difs = ICs - NETs, ptur = ICs - NET0s, pvento = NET0s - NETs)
  e <- df %>% group_by(estacao) %>%
    summarise(across(c(T, UR, V2s, ICs, NETs, difs, ptur, pvento), mean), .groups = "drop")
  tibble(n_est = nrow(e), dif_media = mean(e$difs), parcela_vento = mean(e$pvento),
         parcela_TUR = mean(e$ptur), pct_vento = 100 * mean(e$pvento) / mean(e$difs),
         r_V_dif_estacoes = stats::cor(e$V2s, e$difs), r_V_dif_horarias = stats::cor(df$V2s, df$difs),
         r_pvento_V = stats::cor(e$pvento, e$V2s), r_pvento_T = stats::cor(e$pvento, e$T),
         r_ptur_T = stats::cor(e$ptur, e$T), r_ptur_V = stats::cor(e$ptur, e$V2s))
}
sens <- dplyr::bind_rows(
  decompor(v15) %>% mutate(cenario = "base (V = 0; z0 = 0,03 m)"),
  decompor(v15, excluir_jundiai = TRUE) %>% mutate(cenario = "sem Jundiaí"),
  purrr::map_dfr(P$v_ref_sensibilidade, function(vr) decompor(v15, v_ref = vr) %>%
                   mutate(cenario = paste0("V ref = ", num_br(vr, 1), " m/s"))),
  purrr::map_dfr(P$z0_sensibilidade, function(z) decompor(v15, z0 = z) %>%
                   mutate(cenario = paste0("z0 = ", num_br(z, 2), " m")))) %>% relocate(cenario)
gravar(sens, dG, "sensibilidade_decomposicao.csv")
print(as.data.frame(sens %>% select(cenario, dif_media, parcela_vento, pct_vento)))

# ---- G.6 Jundiaí, temperatura acima de 37 °C, frequência de limiares e derivadas ----------
jund <- val %>% filter(eh_jundiai(estacao)) %>%
  summarise(media_todas_horas = mean(V2), media_15h = mean(V2[hora == P$hora_mapa]),
            maximo = max(V2), pct_menor_0_5 = 100 * mean(V2 < 0.5))
gravar(jund, dG, "jundiai_vento.csv")
acima37 <- tibble(observacoes_15h = nrow(v15), T_acima_37 = sum(v15$T > 37),
                  vento_aquece_NET = sum(v15$p_vento < 0))
gravar(acima37, dG, "observacoes_acima_37C_15h.csv")
freq <- v15 %>% group_by(estacao) %>%
  summarise(n = n(), pct_IC_acima_32 = 100 * mean(IC > P$limiares_frequencia[1]),
            pct_IC_acima_41 = 100 * mean(IC > P$limiares_frequencia[2]),
            ocorrencias_acima_41 = sum(IC > P$limiares_frequencia[2]), .groups = "drop") %>%
  ordenar_estacoes()
gravar(freq, dG, "frequencia_limiares_15h.csv")

# Sensibilidade das fórmulas nas condições médias da rede às 15h (Seção 3.2)
e15 <- est_h[[as.character(P$hora_mapa)]] %>% filter(retida)
Tm <- mean(e15$T); URm <- mean(e15$UR); Vm <- mean(e15$V); hd <- 1e-4
rmse_v <- function(v) met$RMSE[met$variavel == v & met$horario == P$hora_mapa & met$metodo == metodo_vencedor(v, P$hora_mapa)]
derivadas <- tibble(
  T_media = Tm, UR_media = URm, V_media = Vm,
  dIC_dT = (ic_c(Tm + hd, URm) - ic_c(Tm - hd, URm)) / (2 * hd),
  dIC_dUR = (ic_c(Tm, URm + hd) - ic_c(Tm, URm - hd)) / (2 * hd),
  dNET_dT = (net_c(Tm + hd, URm, Vm) - net_c(Tm - hd, URm, Vm)) / (2 * hd),
  dNET_dUR = (net_c(Tm, URm + hd, Vm) - net_c(Tm, URm - hd, Vm)) / (2 * hd)) %>%
  mutate(rmse_T = rmse_v("T"), rmse_UR = rmse_v("UR"),
         erro_IC_por_T = dIC_dT * rmse_T, erro_IC_por_UR = dIC_dUR * rmse_UR,
         erro_NET_por_T = dNET_dT * rmse_T, erro_NET_por_UR = dNET_dUR * rmse_UR)
gravar(derivadas, dG, "derivadas_sensibilidade.csv")

# Gradiente da rede às 15h (Seção 3.1)
gradiente <- with(est_cov[[as.character(P$hora_mapa)]], tibble(
  r_lon_T = stats::cor(lon, T), r_lon_UR = stats::cor(lon, UR), r_distcosta_UR = stats::cor(dist_costa, UR),
  V_min = min(V), V_max = max(V)))
gravar(gradiente, dG, "gradiente_rede_15h.csv")


## ============================================================================
## Etapa H. NET recomposta a partir das variáveis interpoladas
## ============================================================================
msg("Etapa H: NET recomposta")
dH <- pasta("08_recomposta")
recomp <- purrr::map_dfr(P$horas, function(h) {
  pred_de <- function(v) cv %>% filter(horario == h, variavel == v, metodo == metodo_vencedor(v, h)) %>%
    select(estacao, pred)
  d <- cv %>% filter(horario == h, variavel == "NET", metodo == metodo_vencedor("NET", h)) %>%
    select(estacao, NET_obs = obs, NET_direta = pred) %>%
    left_join(pred_de("T") %>% rename(T_pred = pred), by = "estacao") %>%
    left_join(pred_de("UR") %>% rename(UR_pred = pred), by = "estacao") %>%
    left_join(pred_de("V") %>% rename(V_pred = pred), by = "estacao") %>%
    left_join(est_cov[[as.character(h)]] %>% select(estacao, V_era5 = V2m_era5), by = "estacao") %>%
    mutate(UR_pred = pmin(pmax(UR_pred, 0), 100), V_pred = pmax(V_pred, 0),   # limites físicos
           rec_rede = net_c(T_pred, UR_pred, V_pred),
           rec_era5 = net_c(T_pred, UR_pred, V_era5),
           rec_nulo = net_c(T_pred, UR_pred, 0))
  tibble(horario = h, NET_direta = rmse(d$NET_direta - d$NET_obs),
         recomposta_V_rede = rmse(d$rec_rede - d$NET_obs),
         recomposta_V_ERA5 = rmse(d$rec_era5 - d$NET_obs),
         recomposta_V_nulo = rmse(d$rec_nulo - d$NET_obs),
         metodos = paste0("NET ", metodo_vencedor("NET", h), "; T ", metodo_vencedor("T", h), "; UR ",
                          metodo_vencedor("UR", h), "; V ", metodo_vencedor("V", h)))
})
gravar(recomp, dH, "net_recomposta.csv")
print(as.data.frame(recomp))


## ============================================================================
## Etapa I. Efeito do critério de completude (22 estações às 15h)
## ============================================================================
msg("Etapa I: validação cruzada sem o critério de completude")
dI <- pasta("09_completude")
e_todas15 <- est_h[[as.character(P$hora_mapa)]] %>% filter(n > 0)
est22 <- cov_estacoes(pilhas[[as.character(P$hora_mapa)]], e_todas15)
gravar(est22, dI, "estacoes_cov_sem_criterio_15h.csv")
cv22 <- em_cache(file.path(dI, "validacao_cruzada_sem_criterio.rds"), {
  dplyr::bind_rows(validacao_cruzada(est22, "sem critério"), referencias_reanalise(est22)) %>%
    mutate(horario = P$hora_mapa)
}, ass = assinatura(est22))
gravar(cv22, dI, "loocv_sem_criterio_15h.csv")
met22 <- cv22 %>% filter(estacao %in% ret15) %>% metricas()
melhor <- function(m) m %>% filter(metodo %in% METODOS$id, n_falhas == 0) %>% group_by(variavel) %>%
  slice_min(RMSE, n = 1, with_ties = FALSE) %>% ungroup() %>% select(variavel, metodo, RMSE)
contraste <- melhor(met %>% filter(horario == P$hora_mapa)) %>%
  rename(metodo_16 = metodo, RMSE_16 = RMSE) %>%
  left_join(melhor(met22) %>% rename(metodo_22 = metodo, RMSE_22 = RMSE), by = "variavel") %>%
  mutate(variacao_pct = 100 * (RMSE_22 - RMSE_16) / RMSE_16)
v22 <- val %>% filter(hora == P$hora_mapa, estacao %in% e_todas15$estacao)
resumo22 <- tibble(
  conjunto = c("sem critério", "adotado"),
  estacoes = c(nrow(e_todas15), length(ret15)),
  registros_15h = c(nrow(v22), nrow(v15)),
  media_IC = c(mean(e_todas15$IC), mean(e15$IC)), media_NET = c(mean(e_todas15$NET), mean(e15$NET)),
  dif_media = c(mean(e_todas15$IC - e_todas15$NET), mean(e15$IC - e15$NET)))
gravar(contraste, dI, "contraste_16_22_rmse.csv")
gravar(resumo22, dI, "contraste_16_22_resumo.csv")
print(as.data.frame(contraste))


## ============================================================================
## Etapa J. Superfícies de 500 m, distância às estações e zonas de confiança
## ============================================================================
msg("Etapa J: superfícies")
dJ <- pasta("10_superficies")
dJz <- pasta("10_superficies/formato_zenodo")
h_mapa <- as.character(P$hora_mapa)
pilha15 <- pilhas[[h_mapa]]
gdf <- as.data.frame(c(pilha15, mascara), cells = TRUE, na.rm = FALSE)
gdf <- gdf[!is.na(gdf$mascara), ]

superficie <- function(v) {
  cod <- codigo_vencedor(v, P$hora_mapa)
  tr <- dados_variavel(est_cov[[h_mapa]], v)
  nd <- gdf %>% transmute(x, y, altitude, ndvi, lcz, dist_costa, era5 = .data[[COL_ERA5[[v]]]])
  r <- terra::rast(grade); terra::values(r) <- NA_real_
  r[gdf$cell] <- prever(cod, tr, nd)
  names(r) <- v
  r
}
sup <- em_cache(file.path(dJ, "superficies.rds"), {
  s <- terra::rast(purrr::map(VARIAVEIS, superficie))
  terra::wrap(s)
}, ass = assinatura(est_cov[[h_mapa]], venc))
sup <- terra::unwrap(sup)

# Distância de cada célula à estação retida mais próxima (km) e zonas de confiança
pts15 <- terra::vect(pontos_estacoes(e15))
dist_est <- terra::mask(terra::distance(terra::init(grade, 1), pts15) / 1000, mascara)
names(dist_est) <- "dist_estacao_km"
q_dist <- stats::quantile(terra::values(dist_est), c(0.75, 0.90), na.rm = TRUE)
zona <- terra::classify(dist_est, rbind(c(-Inf, q_dist[1], 1), c(q_dist[1], q_dist[2], 2),
                                        c(q_dist[2], Inf, 3)), include.lowest = TRUE)
names(zona) <- "zona_confianca"
dv <- terra::values(dist_est); dv <- dv[!is.na(dv)]
area_km2 <- as.numeric(sf::st_area(rn)) / 1e6
resumo_dist <- tibble(
  area_km2 = area_km2, km2_por_estacao = area_km2 / nrow(e15),
  dist_media_km = mean(dv), p75_km = unname(q_dist[1]), p90_km = unname(q_dist[2]),
  dist_max_km = max(dv), pct_area_acima_50km = 100 * mean(dv > 50))
gravar(resumo_dist, dJ, "distancia_estacoes_resumo.csv")

# Preenchimento só para exibição: células sem covariável dentro do limite sem lagoas
# e salinas recebem a média das vizinhas (Seção 3.8).
preencher <- function(r) {
  alvo <- !is.na(terra::values(mascara_exib))
  for (k in 1:60) {
    if (!any(is.na(terra::values(r)[alvo]))) break
    r <- terra::focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
  }
  terra::mask(r, mascara_exib)
}
sup_exib <- terra::rast(lapply(seq_len(terra::nlyr(sup)), function(i) preencher(sup[[i]])))
names(sup_exib) <- names(sup)

for (v in VARIAVEIS) {
  terra::writeRaster(sup[[v]], file.path(dJ, sprintf("RN_%s_15h_2025_analise.tif", v)),
                     overwrite = TRUE, datatype = "FLT4S")
  terra::writeRaster(sup_exib[[v]], file.path(dJz, sprintf("RN_%s_15h_2025.tif", v)),
                     overwrite = TRUE, datatype = "FLT4S")
}
terra::writeRaster(dist_est, file.path(dJz, "RN_dist_estacao_km_2025.tif"), overwrite = TRUE, datatype = "FLT4S")
terra::writeRaster(zona, file.path(dJz, "RN_zona_confianca_15h_2025.tif"), overwrite = TRUE, datatype = "INT1U")

faixas_grade <- purrr::map_dfr(VARIAVEIS, function(v) {
  x <- terra::values(sup_exib[[v]]); x <- x[!is.na(x)]
  tibble(variavel = v, metodo = metodo_vencedor(v, P$hora_mapa), minimo = min(x), maximo = max(x),
         celulas = length(x),
         minimo_estacoes = min(e15[[v]]), maximo_estacoes = max(e15[[v]]))
})
gravar(faixas_grade, dJ, "faixas_das_superficies.csv")

# Estações com as predições da validação cruzada dos vencedores (formato do Zenodo)
est_zen <- e15 %>% select(estacao, rede, lon, lat, n_registros = n, T, UR, V, IC, NET) %>%
  left_join(est_cov[[h_mapa]] %>% select(estacao, altitude, dist_costa_km = dist_costa), by = "estacao")
for (v in VARIAVEIS) {
  est_zen <- est_zen %>% left_join(cv %>% filter(horario == P$hora_mapa, variavel == v,
                                                 metodo == metodo_vencedor(v, P$hora_mapa)) %>%
                                     select(estacao, !!paste0(v, "_loocv") := pred), by = "estacao")
}
gravar(est_zen, dJz, "estacoes_15h_2025.csv")
gravar(met %>% filter(horario == P$hora_mapa), dJz, "validacao_metodos_15h_2025.csv")


## ============================================================================
## Etapa K. Leitura municipal, população e escalas de classificação
## ============================================================================
msg("Etapa K: leitura municipal")
dK <- pasta("11_municipios")
camadas <- c(sup, dist_est)
ex <- terra::extract(camadas, terra::vect(mun), exact = TRUE)
mm <- ex %>% group_by(ID) %>%
  summarise(across(all_of(names(camadas)), ~ sum(.x * fraction, na.rm = TRUE) / sum(fraction[!is.na(.x)])),
            .groups = "drop")
munt <- sf::st_drop_geometry(mun) %>% mutate(ID = dplyr::row_number()) %>% left_join(mm, by = "ID")

if (existe(ARQ_POPULACAO)) {
  pop <- readr::read_csv(ARQ_POPULACAO, show_col_types = FALSE)
  col_cod <- c(COLS_POP[["codigo"]], "code_muni", "cod_ibge", "CD_MUN", "codigo", "cod_mun")
  col_pop <- c(COLS_POP[["populacao"]], "populacao", "pop", "pop_2022", "populacao_2022", "Populacao", "POP")
  col_cod <- intersect(col_cod, names(pop))[1]; col_pop <- intersect(col_pop, names(pop))[1]
  if (is.na(col_cod) || is.na(col_pop)) {
    stop("Colunas de código e população não encontradas em ", basename(ARQ_POPULACAO), ". Disponíveis: ",
         paste(names(pop), collapse = ", "), "\nAjuste COLS_POP.", call. = FALSE)
  }
  pop <- tibble(code_muni = num(pop[[col_cod]]), populacao = num(pop[[col_pop]])) %>%
    distinct(code_muni, .keep_all = TRUE)
  munt <- munt %>% left_join(pop, by = "code_muni")
  if (anyNA(munt$populacao)) warning(sum(is.na(munt$populacao)), " municípios sem população.", call. = FALSE)
} else {
  warning("Arquivo de população não encontrado; as tabelas populacionais ficarão vazias.", call. = FALSE)
  munt$populacao <- NA_real_
}

rmse_ic <- met$RMSE[met$variavel == "IC" & met$horario == P$hora_mapa & met$metodo == metodo_vencedor("IC", P$hora_mapa)]
rot_nws <- c("Seguro (até 27 °C)", "Cuidado (27 a 32 °C)", "Cuidado extremo (32 a 41 °C)",
             "Perigo (41 a 54 °C)", "Perigo extremo (acima de 54 °C)")
rot_for <- c("Não designado (abaixo de 29,3 °C)", "Confortável (29,3 a 32,1 °C)",
             "Desconfortável (acima de 32,1 °C)")
munt <- munt %>% mutate(
  dif_IC_NET = IC - NET,
  classe_NWS = cut(IC, c(-Inf, P$limiares_nws, Inf), labels = rot_nws, right = TRUE),
  classe_Fortaleza = cut(IC, c(-Inf, P$limiares_fortaleza, Inf), labels = rot_for, right = TRUE),
  dist_limiar_NWS = pmin(abs(IC - P$limiares_nws[1]), abs(IC - P$limiares_nws[2]), abs(IC - P$limiares_nws[3])),
  perto_de_limiar = dist_limiar_NWS < rmse_ic,
  fora_zona_alta = dist_estacao_km > q_dist[1],
  posicao_IC = rank(-IC, ties.method = "first"), posicao_NET = rank(-NET, ties.method = "first"),
  deslocamento = posicao_IC - posicao_NET)
gravar(munt, dK, "municipios_15h_2025.csv")
gravar(munt %>% select(code_muni, name_muni, IC, NET, T, UR, V, dist_estacao_km, populacao,
                       classe_NWS, classe_Fortaleza, posicao_IC, posicao_NET, deslocamento),
       dJz, "municipios_15h_2025.csv")

pop_tot <- sum(munt$populacao, na.rm = TRUE)
classif <- dplyr::bind_rows(
  munt %>% count(classe = classe_NWS, wt = NULL, name = "municipios") %>%
    left_join(munt %>% group_by(classe = classe_NWS) %>% summarise(populacao = sum(populacao)), by = "classe") %>%
    mutate(escala = "National Weather Service"),
  munt %>% count(classe = classe_Fortaleza, name = "municipios") %>%
    left_join(munt %>% group_by(classe = classe_Fortaleza) %>% summarise(populacao = sum(populacao)), by = "classe") %>%
    mutate(escala = "Ajustada para Fortaleza")) %>%
  mutate(classe = as.character(classe), pct_populacao = 100 * populacao / pop_tot) %>%
  relocate(escala)
gravar(classif, dK, "classificacao_municipal.csv")
print(as.data.frame(classif))

leitura <- tibble(
  municipios = nrow(munt), populacao_total = pop_tot,
  r_IC_NET_municipal = stats::cor(munt$IC, munt$NET),
  dif_min = min(munt$dif_IC_NET), dif_max = max(munt$dif_IC_NET),
  desloc_mais_30 = sum(abs(munt$deslocamento) > 30), desloc_mais_20 = sum(abs(munt$deslocamento) > 20),
  r_deslocamento_V = stats::cor(munt$deslocamento, munt$V),
  IC_mun_min = min(munt$IC), IC_mun_max = max(munt$IC),
  NET_mun_fora_faixa_obs = sum(munt$NET < min(e15$NET) | munt$NET > max(e15$NET)),
  rmse_IC_vencedor = rmse_ic,
  perto_limiar_municipios = sum(munt$perto_de_limiar),
  perto_limiar_pct_pop = 100 * sum(munt$populacao[munt$perto_de_limiar]) / pop_tot,
  fora_zona_alta_municipios = sum(munt$fora_zona_alta),
  fora_zona_alta_populacao = sum(munt$populacao[munt$fora_zona_alta]),
  municipios_acima_40km = sum(munt$dist_estacao_km > 40),
  diferem_NWS_extremo_vs_Fortaleza_desconf = sum(xor(munt$IC > P$limiares_nws[2], munt$IC > P$limiares_fortaleza[2])))
gravar(leitura, dK, "leitura_municipal_resumo.csv")


## ============================================================================
## Etapa L. Estações de Natal (INMET A304 e PCD Natal)
## ============================================================================
msg("Etapa L: estações de Natal")
dL <- pasta("12_natal")
nat <- val %>% filter(grepl("A304", estacao) | (rede == "EMPARN" & chave(estacao) == "PCD NATAL"))
if (length(unique(nat$rede)) == 2) {
  nat15 <- nat %>% filter(hora == P$hora_mapa) %>%
    select(rede, data, T, UR, IC, V_nativo = V, V2) %>%
    pivot_wider(names_from = rede, values_from = c(T, UR, IC, V_nativo, V2)) %>% tidyr::drop_na()
  natal <- tibble(
    variavel = c("T (°C)", "UR (%)", "IC (°C)", "V na altura nativa (m/s)", "V a 2 m (m/s)"),
    INMET = c(mean(nat15$T_INMET), mean(nat15$UR_INMET), mean(nat15$IC_INMET),
              mean(nat15$V_nativo_INMET), mean(nat15$V2_INMET)),
    EMPARN = c(mean(nat15$T_EMPARN), mean(nat15$UR_EMPARN), mean(nat15$IC_EMPARN),
               mean(nat15$V_nativo_EMPARN), mean(nat15$V2_EMPARN))) %>%
    mutate(diferenca_INMET_menos_EMPARN = INMET - EMPARN, dias_comuns = nrow(nat15))
  natal_hora <- nat %>% select(rede, datahora_utc, hora, T) %>%
    pivot_wider(names_from = rede, values_from = T) %>% tidyr::drop_na() %>%
    group_by(hora) %>% summarise(n = n(), dif_T_INMET_menos_EMPARN = mean(INMET - EMPARN), .groups = "drop")
  gravar(natal, dL, "natal_dias_comuns_15h.csv")
  gravar(natal_hora, dL, "natal_diferenca_T_por_hora.csv")
  print(as.data.frame(natal))
} else {
  natal <- NULL
  msg("   as duas estações de Natal não foram encontradas; etapa ignorada")
}


## ============================================================================
## Etapa M. Tabelas do artigo (CSV com ; e vírgula decimal, para o Word ou Excel)
## ============================================================================
msg("Etapa M: tabelas")
dM <- pasta("13_tabelas")

# Tabela 1: estações retidas às 15h
t1 <- e15 %>% ordenar_estacoes() %>%
  transmute(`Estação` = rotulo_estacao(estacao), Rede = rede, n = n,
            `T (°C)` = T, `UR (%)` = UR, `V (m/s)` = V, `IC (°C)` = IC, `NET (°C)` = NET,
            `IC - NET (°C)` = IC - NET)
t1 <- dplyr::bind_rows(t1, t1 %>% summarise(`Estação` = "Média", Rede = "", n = mean(n),
                                            across(`T (°C)`:`IC - NET (°C)`, mean)))
t1f <- t1 %>% mutate(n = num_br(n, 0), across(c(`T (°C)`, `UR (%)`, `IC (°C)`, `NET (°C)`, `IC - NET (°C)`), ~ num_br(.x, 1)),
                     `V (m/s)` = num_br(`V (m/s)`, 2))
gravar_br(t1f, dM, "Tabela1_estacoes_15h.csv")

# Tabela 2: R² e RMSE às 15h por método e variável
t2 <- met %>% filter(horario == P$hora_mapa) %>%
  mutate(ordem = match(metodo, ROTULOS$id), R2f = fmt_r2(R2), RMSEf = fmt_rmse(RMSE, variavel)) %>%
  select(ordem, metodo, rotulo, variavel, R2f, RMSEf) %>%
  pivot_wider(names_from = variavel, values_from = c(R2f, RMSEf)) %>% arrange(ordem)
t2f <- t2 %>% transmute(`Método` = paste(metodo, rotulo),
                        `IC R²` = R2f_IC, `IC RMSE` = RMSEf_IC, `NET R²` = R2f_NET, `NET RMSE` = RMSEf_NET,
                        `T R²` = R2f_T, `T RMSE` = RMSEf_T, `UR R²` = R2f_UR, `UR RMSE` = RMSEf_UR,
                        `V R²` = R2f_V, `V RMSE` = RMSEf_V)
gravar_br(t2f, dM, "Tabela2_validacao_15h.csv")

# Tabela 3: vencedores nos três horários
t3 <- venc %>% mutate(celula = paste0(ifelse(sem_poder_preditivo, "sem poder preditivo: ", ""),
                                      metodo, " (", fmt_r2(R2), "; ", fmt_rmse(RMSE, variavel), ")"),
                      horario = sprintf("%02dh", horario)) %>%
  select(variavel, horario, celula) %>% pivot_wider(names_from = horario, values_from = celula) %>%
  mutate(variavel = factor(variavel, VARIAVEIS)) %>% arrange(variavel) %>% rename(`Variável` = variavel)
gravar_br(t3, dM, "Tabela3_vencedores.csv")

# Quadro 1: métodos avaliados
q1 <- ROTULOS %>% transmute(ID = id, `Família` = familia, `Método` = rotulo, `Código` = codigo,
                            `Covariáveis` = vapply(codigo, function(cc) {
                              v <- NECESSARIAS[[cc]]
                              if (is.null(v)) "ERA5" else paste(v, collapse = ", ")
                            }, character(1)))
gravar_br(q1, dM, "Quadro1_metodos.csv")

# Tabela 4: classes de vento às 15h
t4 <- classes %>% transmute(`V (m/s)` = as.character(classe), n = num_br(n, 0), `%` = num_br(pct, 1),
                            `T (°C)` = num_br(T, 1), `UR (%)` = num_br(UR, 1), `IC (°C)` = num_br(IC, 1),
                            `NET (°C)` = num_br(NET, 1), `IC - NET` = num_br(dif, 2),
                            `Parcela T e UR` = num_br(p_tur, 2), `Parcela do vento` = num_br(p_vento, 2))
gravar_br(t4, dM, "Tabela4_classes_vento_15h.csv")

# Tabela 5: NET recomposta
pct_txt <- function(a, b) paste0(num_br(a, 3), " (", ifelse(a >= b, "+", "\u2212"),
                                 num_br(abs(100 * (a - b) / b), 1), "%)")
t5 <- recomp %>% transmute(`Horário` = sprintf("%02dh", horario), `NET direta` = num_br(NET_direta, 3),
                           `Recomposta, V da rede` = pct_txt(recomposta_V_rede, NET_direta),
                           `Recomposta, V do ERA5-Land` = pct_txt(recomposta_V_ERA5, NET_direta),
                           `Recomposta, V nulo` = pct_txt(recomposta_V_nulo, NET_direta))
gravar_br(t5, dM, "Tabela5_NET_recomposta.csv")

# Tabela 6: classificação municipal
t6 <- classif %>% filter(municipios > 0) %>%
  transmute(Escala = escala, Faixa = classe, `Municípios` = num_br(municipios, 0),
            `População` = num_br(populacao, 0), `% da população` = num_br(pct_populacao, 1))
gravar_br(t6, dM, "Tabela6_classificacao_municipal.csv")

# Tabelas S1 e S2 (material suplementar)
tS1 <- kg$parametros %>% mutate(`Estação` = rotulo_estacao(estacao)) %>%
  transmute(Variante = variante, `Estação`, Modelo = modelo, Pepita = num_br(pepita, 3),
            `Patamar parcial` = num_br(patamar_parcial, 3), `Alcance (km)` = num_br(alcance_m / 1000, 1),
            Observado = num_br(obs, 2), Predito = num_br(pred, 2))
gravar_br(tS1, dM, "TabelaS1_krigagem.csv")
tS2 <- inc$ic_rmse %>% mutate(ordem = match(metodo, ROTULOS$id)) %>%
  arrange(factor(variavel, VARIAVEIS), horario, ordem) %>%
  transmute(`Variável` = variavel, `Horário` = sprintf("%02dh", horario), `Método` = metodo,
            n, RMSE = fmt_rmse(rmse, variavel), `IC95% inferior` = fmt_rmse(ic95_inf, variavel),
            `IC95% superior` = fmt_rmse(ic95_sup, variavel))
gravar_br(tS2, dM, "TabelaS2_rmse_ic95.csv")


## ============================================================================
## Etapa N. Figuras do artigo (legíveis em tons de cinza)
## ============================================================================
if (GERAR_FIGURAS) {
  msg("Etapa N: figuras")
  dN <- pasta("14_figuras")
  tema <- theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.tag = element_text(face = "bold"))
  salvar_fig <- function(p, nome, w = 16, h = 10) {
    ggsave(file.path(dN, nome), p, width = w, height = h, units = "cm", dpi = 300, bg = "white")
  }
  juntar <- function(a, b, ncol = 2) {
    if (tem("patchwork")) return(patchwork::wrap_plots(a, b, ncol = ncol) +
                                   patchwork::plot_annotation(tag_levels = list(c("(a)", "(b)"))))
    NULL
  }
  coord_mapa <- coord_sf(crs = CRS_GRADE, datum = sf::st_crs(4326), expand = FALSE)
  rn_sf <- sf::st_sf(geometry = rn_cheio)
  mun_sf <- mun

  # Figura 1: rede e critério de completude às 15h
  st1 <- pontos_estacoes(est_h[[h_mapa]]) %>% filter(!is.na(lon), !is.na(lat)) %>%
    mutate(rotulo = rotulo_estacao(estacao, "figura"),
           situacao = factor(situacao, c("Retida na análise", paste0("Abaixo de ", P$min_registros, " registros"),
                                         "Sem registro conjunto das três variáveis")))
  cont <- table(st1$situacao)
  levels(st1$situacao) <- paste0(levels(st1$situacao), " (", cont[levels(st1$situacao)], ")")
  xy1 <- sf::st_coordinates(st1)
  st1$X <- xy1[, 1]; st1$Y <- xy1[, 2]
  f1 <- ggplot() +
    geom_sf(data = mun_sf, fill = "grey96", colour = "grey80", linewidth = 0.15) +
    geom_sf(data = rn_sf, fill = NA, colour = "grey25", linewidth = 0.4) +
    geom_point(data = st1, aes(X, Y, shape = situacao, fill = situacao), size = 2.6, stroke = 0.7) +
    scale_shape_manual(values = c(21, 24, 4), name = NULL) +
    scale_fill_manual(values = c("black", "white", "black"), name = NULL) +
    coord_sf(crs = CRS_GRADE, datum = sf::st_crs(4326), expand = TRUE, clip = "off") +
    tema + theme(legend.direction = "vertical") + labs(x = NULL, y = NULL)
  f1 <- if (tem("ggrepel")) {
    f1 + ggrepel::geom_text_repel(data = st1, aes(X, Y, label = rotulo), size = 2.4, min.segment.length = 0.2,
                                  seed = 1, max.overlaps = Inf)
  } else f1 + geom_text(data = st1, aes(X, Y, label = rotulo), size = 2.3, vjust = -0.9)
  salvar_fig(f1, "Figura1_rede_estacoes.png", 16, 11)

  # Figura 2: RMSE da validação cruzada às 15h (escala logarítmica)
  d2 <- met %>% filter(horario == P$hora_mapa, metodo %in% METODOS$id) %>%
    mutate(rotulo_m = factor(paste(metodo, rotulo), rev(paste(METODOS$id, METODOS$rotulo))),
           vencedor = paste(variavel, metodo) %in% paste(venc$variavel[venc$horario == P$hora_mapa],
                                                          venc$metodo[venc$horario == P$hora_mapa]),
           painel = factor(paste0(variavel, " (", UNIDADE[variavel], ")"), paste0(VARIAVEIS, " (", UNIDADE, ")")))
  r2_lin <- met %>% filter(horario == P$hora_mapa, metodo == "R2") %>%
    mutate(painel = factor(paste0(variavel, " (", UNIDADE[variavel], ")"), levels(d2$painel)))
  f2 <- ggplot(d2, aes(RMSE, rotulo_m)) +
    geom_vline(data = r2_lin, aes(xintercept = RMSE), linetype = "dashed", colour = "grey30") +
    geom_point(aes(shape = vencedor, fill = vencedor), size = 2.2) +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
    scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "black"), guide = "none") +
    scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000),
                  labels = function(x) ifelse(x < 1, num_br(x, 1), num_br(x, 0))) +
    facet_wrap(~ painel, nrow = 1, scales = "free_x") +
    labs(x = "RMSE da validação cruzada (escala logarítmica)", y = NULL) + tema +
    theme(panel.spacing.x = unit(0.9, "lines"), axis.text.x = element_text(size = 7))
  salvar_fig(f2, "Figura2_rmse_metodos_15h.png", 19, 10)

  # Figura 3: decomposição às 15h
  long_cl <- classes %>% select(V_media, dif, p_tur, p_vento) %>%
    pivot_longer(-V_media, names_to = "serie", values_to = "valor") %>%
    mutate(serie = factor(serie, c("dif", "p_tur", "p_vento"),
                          c("Diferença IC - NET", "Parcela de T e UR", "Parcela do vento")))
  f3a <- ggplot() +
    { if (tem("hexbin")) geom_hex(data = v15, aes(V2, dif), bins = 35, fill = "grey85", colour = NA)
      else geom_bin2d(data = v15, aes(V2, dif), bins = 35, fill = "grey85") } +
    geom_line(data = long_cl, aes(V_media, valor, linetype = serie), colour = "black") +
    geom_point(data = long_cl, aes(V_media, valor, shape = serie), size = 2, fill = "white") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = NULL) +
    scale_shape_manual(values = c(16, 21, 17), name = NULL) +
    labs(x = "Velocidade do vento a 2 m (m/s)", y = "°C") + tema + theme(legend.direction = "vertical")
  d3b <- dec_est %>% arrange(desc(dist_costa)) %>%
    mutate(rot = paste0(rotulo_estacao(estacao, "figura"), " [", num_br(V2, 1), "]"),
           rot = factor(rot, rot)) %>%
    select(rot, p_vento, p_tur) %>% pivot_longer(-rot, names_to = "parcela", values_to = "valor") %>%
    mutate(parcela = factor(parcela, c("p_tur", "p_vento"), c("Parcela de T e UR", "Parcela do vento")))
  f3b <- ggplot(d3b, aes(valor, rot, fill = parcela)) + geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
    scale_fill_manual(values = c("grey80", "grey25"), name = NULL) +
    labs(x = "Diferença IC - NET (°C)", y = NULL) + tema
  f3 <- juntar(f3a, f3b)
  if (!is.null(f3)) salvar_fig(f3, "Figura3_decomposicao_15h.png", 19, 10) else {
    salvar_fig(f3a, "Figura3a_decomposicao_15h.png", 10, 10); salvar_fig(f3b, "Figura3b_decomposicao_15h.png", 10, 10)
  }

  # Figura 4: ciclo diário
  l4 <- ciclo %>% select(hora, dif, p_tur, p_vento) %>% pivot_longer(-hora, names_to = "serie", values_to = "valor") %>%
    mutate(serie = factor(serie, c("dif", "p_tur", "p_vento"),
                          c("Diferença IC - NET", "Parcela de T e UR", "Parcela do vento")))
  f4a <- ggplot(l4, aes(hora, valor, linetype = serie, shape = serie)) + geom_line() +
    geom_point(size = 1.6, fill = "white") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = NULL) +
    scale_shape_manual(values = c(16, 21, 17), name = NULL) +
    scale_x_continuous(breaks = seq(0, 21, 3)) + labs(x = "Hora local", y = "°C") + tema +
    theme(legend.direction = "vertical")
  f4b <- ggplot(ciclo, aes(hora, V2)) + geom_line() + geom_point(size = 1.6) +
    scale_x_continuous(breaks = seq(0, 21, 3)) + labs(x = "Hora local", y = "V a 2 m (m/s)") + tema
  f4 <- juntar(f4a, f4b)
  if (!is.null(f4)) salvar_fig(f4, "Figura4_ciclo_diario.png", 18, 8.5) else {
    salvar_fig(f4a, "Figura4a_ciclo_diario.png", 10, 8.5); salvar_fig(f4b, "Figura4b_vento.png", 9, 8.5)
  }

  # Figura 5: superfícies de IC e NET às 15h
  mapa_raster <- function(r, titulo) {
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(df)[3] <- "v"
    ggplot() + geom_raster(data = df, aes(x, y, fill = v)) +
      geom_sf(data = mun_sf, fill = NA, colour = "white", linewidth = 0.1) +
      geom_sf(data = rn_sf, fill = NA, colour = "grey20", linewidth = 0.3) +
      geom_sf(data = pontos_estacoes(e15), shape = 21, fill = "white", size = 1.4) +
      scale_fill_distiller(palette = "YlOrRd", direction = 1, name = titulo) +
      coord_mapa + labs(x = NULL, y = NULL) + tema +
      theme(legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.25, "cm"))
  }
  f5a <- mapa_raster(sup_exib[["IC"]], "IC (°C)")
  f5b <- mapa_raster(sup_exib[["NET"]], "NET (°C)")
  f5 <- juntar(f5a, f5b, ncol = 1)
  if (!is.null(f5)) salvar_fig(f5, "Figura5_superficies_15h.png", 12, 17) else {
    salvar_fig(f5a, "Figura5a_IC_15h.png", 12, 8.5); salvar_fig(f5b, "Figura5b_NET_15h.png", 12, 8.5)
  }

  # Figura 6: leitura municipal (diferença e deslocamento de posição)
  mun_map <- mun_sf %>% left_join(munt %>% select(code_muni, dif_IC_NET, deslocamento), by = "code_muni")
  f6a <- ggplot(mun_map) + geom_sf(aes(fill = dif_IC_NET), colour = "white", linewidth = 0.1) +
    geom_sf(data = rn_sf, fill = NA, colour = "grey20", linewidth = 0.3) +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "IC - NET (°C)") + coord_mapa + tema +
    theme(legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.25, "cm"))
  lim6 <- max(abs(mun_map$deslocamento), na.rm = TRUE)
  f6b <- ggplot(mun_map) + geom_sf(aes(fill = deslocamento), colour = "white", linewidth = 0.1) +
    geom_sf(data = rn_sf, fill = NA, colour = "grey20", linewidth = 0.3) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-lim6, lim6),
                         name = "Deslocamento (posições)") + coord_mapa + tema +
    theme(legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.25, "cm"))
  f6 <- juntar(f6a, f6b, ncol = 1)
  if (!is.null(f6)) salvar_fig(f6, "Figura6_leitura_municipal.png", 12, 17) else {
    salvar_fig(f6a, "Figura6a_diferenca.png", 12, 8.5); salvar_fig(f6b, "Figura6b_deslocamento.png", 12, 8.5)
  }

  # Figura 7: classificação municipal do IC pelas duas escalas
  mun_cl <- mun_sf %>% left_join(munt %>% select(code_muni, classe_NWS, classe_Fortaleza, fora_zona_alta),
                                 by = "code_muni")
  fora <- mun_cl %>% filter(fora_zona_alta)
  cores_classe <- c(stats::setNames(c("grey97", "grey80", "grey50", "grey25", "black"), rot_nws),
                    stats::setNames(c("grey97", "grey80", "grey50"), rot_for))
  mapa_classe <- function(campo) {
    ggplot(mun_cl) + geom_sf(aes(fill = droplevels(.data[[campo]])), colour = "white", linewidth = 0.1) +
      geom_sf(data = fora, fill = NA, colour = "black", linewidth = 0.45) +
      geom_sf(data = rn_sf, fill = NA, colour = "grey20", linewidth = 0.3) +
      scale_fill_manual(values = cores_classe, name = NULL, drop = TRUE) +
      coord_mapa + tema + theme(legend.direction = "vertical")
  }
  f7 <- juntar(mapa_classe("classe_NWS"), mapa_classe("classe_Fortaleza"), ncol = 1)
  if (!is.null(f7)) salvar_fig(f7, "Figura7_classificacao_municipal.png", 12, 17) else {
    salvar_fig(mapa_classe("classe_NWS"), "Figura7a_NWS.png", 12, 9)
    salvar_fig(mapa_classe("classe_Fortaleza"), "Figura7b_Fortaleza.png", 12, 9)
  }
}


## ============================================================================
## Etapa O. Conferência com os valores publicados no artigo
## ============================================================================
msg("Etapa O: conferência com os valores publicados")
dO <- pasta("15_conferencia")
m15 <- function(v, m, col = "RMSE") {
  x <- met[[col]][met$variavel == v & met$horario == P$hora_mapa & met$metodo == m]
  if (length(x)) x else NA_real_
}
r_h <- function(h, col) { x <- resumo_rede[[col]][resumo_rede$horario == h]; if (length(x)) x else NA_real_ }
kgv <- function(v, col) { x <- kg_resumo[[col]][kg_resumo$variante == v]; if (length(x)) x else NA_real_ }
cls <- function(esc, padrao, col) {
  x <- classif[[col]][classif$escala == esc & grepl(padrao, classif$classe)]
  if (length(x)) sum(x) else 0
}
rec15 <- function(col) recomp[[col]][recomp$horario == P$hora_mapa]
fx <- function(v, col) faixas_grade[[col]][faixas_grade$variavel == v]
kauto <- kg$parametros %>% filter(variante == "a_auto_por_rodada")
L <- function(item, publicado, tolerancia, obtido) {
  tibble(item = item, publicado = publicado, tolerancia = tolerancia,
         obtido = if (length(obtido) == 1) as.numeric(obtido) else NA_real_)
}
conf <- dplyr::bind_rows(
  L("Células da grade dentro do limite do RN (análise)", 212232, 0, info_grade$celulas_analise),
  L("Células da grade com valor nas superfícies publicadas", 213600, 0, info_grade$celulas_exibicao),
  L("Estações retidas às 15h", 16, 0, r_h(15, "estacoes_retidas")),
  L("Estações retidas às 09h", 16, 0, r_h(9, "estacoes_retidas")),
  L("Estações retidas às 00h", 14, 0, r_h(0, "estacoes_retidas")),
  L("Estações abaixo do limiar às 15h", 6, 0, r_h(15, "abaixo_do_limiar")),
  L("Estações sem registro conjunto às 15h", 3, 0, r_h(15, "sem_registro")),
  L("Registros válidos às 15h (estações retidas)", 5201, 0, r_h(15, "registros_retidas")),
  L("Registros válidos às 09h (estações retidas)", 4923, 0, r_h(9, "registros_retidas")),
  L("Registros válidos às 00h (estações retidas)", 4403, 0, r_h(0, "registros_retidas")),
  L("Média do IC nas estações às 15h (°C)", 32.8, 0.05, mean(e15$IC)),
  L("Média da NET nas estações às 15h (°C)", 25.1, 0.05, mean(e15$NET)),
  L("r(longitude, T) às 15h", -0.85, 0.005, gradiente$r_lon_T),
  L("r(longitude, UR) às 15h", 0.80, 0.005, gradiente$r_lon_UR),
  L("r(distância à costa, UR) às 15h", -0.87, 0.005, gradiente$r_distcosta_UR),
  L("Diferença média IC - NET às 15h (°C)", 7.697, 0.0005, sens$dif_media[1]),
  L("Parcela do vento às 15h (°C)", 2.370, 0.0005, sens$parcela_vento[1]),
  L("Parcela do vento às 15h (%)", 30.79, 0.005, sens$pct_vento[1]),
  L("r(V, IC - NET) entre estações", 0.777, 0.0005, sens$r_V_dif_estacoes[1]),
  L("r(V, IC - NET) entre observações das 15h", 0.534, 0.0005, sens$r_V_dif_horarias[1]),
  L("r(parcela do vento, V) entre estações", 0.73, 0.005, sens$r_pvento_V[1]),
  L("r(parcela do vento, T) entre estações", -0.79, 0.005, sens$r_pvento_T[1]),
  L("r(parcela de T e UR, T) entre estações", 0.87, 0.005, sens$r_ptur_T[1]),
  L("Sem Jundiaí: diferença média (°C)", 7.93, 0.005, sens$dif_media[2]),
  L("Sem Jundiaí: parcela do vento (%)", 30.85, 0.005, sens$pct_vento[2]),
  L("Sem Jundiaí: r(V, IC - NET) entre estações", 0.61, 0.005, sens$r_V_dif_estacoes[2]),
  L("V ref 0,1 m/s: parcela do vento (%)", 26.30, 0.005, sens$pct_vento[3]),
  L("V ref 0,2 m/s: parcela do vento (%)", 23.56, 0.005, sens$pct_vento[4]),
  L("z0 = 0,01 m: parcela do vento (%)", 30.93, 0.005, sens$pct_vento[5]),
  L("z0 = 0,25 m: parcela do vento (%)", 30.24, 0.005, sens$pct_vento[7]),
  L("Observações das 15h com T acima de 37 °C", 78, 0, acima37$T_acima_37),
  L("Observações de todas as horas (todas as estações)", 126698, 0,
    todas$observacoes[todas$conjunto == "todas as estações"]),
  L("r(V, IC - NET) em todas as horas", 0.75, 0.005, todas$r_V_dif[todas$conjunto == "todas as estações"]),
  L("Observações de calmaria (todas as horas, todas as estações)", 9121, 0,
    calma$observacoes[calma$conjunto == "todas as estações, todas as horas"]),
  L("Observações de calmaria às 15h", 215, 0, calma$observacoes[calma$conjunto == "16 estações retidas, 15h"]),
  L("Jundiaí: vento médio em todas as horas (m/s)", 0.172, 0.0005, if (nrow(jund)) jund$media_todas_horas else NA),
  L("Jundiaí: vento máximo (m/s)", 1.04, 0.005, if (nrow(jund)) jund$maximo else NA),
  L("RMSE M10, IC, 15h", 1.471, 0.0005, m15("IC", "M10")),
  L("R² M10, IC, 15h", 0.529, 0.0005, m15("IC", "M10", "R2")),
  L("RMSE M11, NET, 15h", 0.855, 0.0005, m15("NET", "M11")),
  L("R² M11, NET, 15h", 0.758, 0.0005, m15("NET", "M11", "R2")),
  L("RMSE M10, T, 15h", 0.761, 0.0005, m15("T", "M10")),
  L("RMSE M14, UR, 15h", 4.26, 0.005, m15("UR", "M14")),
  L("RMSE M10, V, 15h", 1.000, 0.0005, m15("V", "M10")),
  L("RMSE M1 (IDW), IC, 15h", 1.540, 0.0005, m15("IC", "M1")),
  L("RMSE M5 (coordenadas), IC, 15h", 1.814, 0.0005, m15("IC", "M5")),
  L("RMSE M6 (coordenadas e altitude), IC, 15h", 1.769, 0.0005, m15("IC", "M6")),
  L("RMSE M8 (GAM coordenadas), IC, 15h", 1.918, 0.0005, m15("IC", "M8")),
  L("RMSE M13 (floresta aleatória), IC, 15h", 1.563, 0.0005, m15("IC", "M13")),
  L("RMSE M2 (krigagem ordinária), IC, 15h", 198, 0.5, m15("IC", "M2")),
  L("RMSE R1, IC, 15h", 1.450, 0.0005, m15("IC", "R1")),
  L("RMSE R2, IC, 15h", 1.346, 0.0005, m15("IC", "R2")),
  L("RMSE R2, NET, 15h", 0.784, 0.0005, m15("NET", "R2")),
  L("RMSE R2, UR, 15h", 4.20, 0.005, m15("UR", "R2")),
  L("RMSE R2, V, 15h", 0.943, 0.0005, m15("V", "R2")),
  L("Krigagem automática: rodadas com modelo de Matérn (Ste)", 15, 0, sum(kauto$modelo == "Ste", na.rm = TRUE)),
  L("Krigagem automática: rodadas com pepita nula", 14, 0, kgv("a_auto_por_rodada", "pepita_nula")),
  L("Krigagem automática: alcance mediano (km)", 31.8, 0.05, kgv("a_auto_por_rodada", "alcance_mediano_km")),
  L("Krigagem automática: predições a mais de 1 °C da faixa", 11, 0, kgv("a_auto_por_rodada", "fora_da_faixa_mais_1C")),
  L("Krigagem com variograma fixo: RMSE do IC", 1.567, 0.0005, kgv("b_fixo_todas", "rmse")),
  L("Krigagem com variograma fixo: alcance (km)", 780, 0.5, kgv("b_fixo_todas", "alcance_mediano_km")),
  L("Krigagem por REML: estações preditas", 2, 0, kgv("d_REML_todas", "n_predicoes")),
  L("NET recomposta às 15h, V da rede (RMSE)", 1.040, 0.0005, rec15("recomposta_V_rede")),
  L("NET recomposta às 15h, V do ERA5-Land (RMSE)", 1.060, 0.0005, rec15("recomposta_V_ERA5")),
  L("NET recomposta às 15h, V nulo (RMSE)", 2.784, 0.0005, rec15("recomposta_V_nulo")),
  L("Sem critério: estações às 15h", 22, 0, resumo22$estacoes[1]),
  L("Sem critério: registros às 15h", 5838, 0, resumo22$registros_15h[1]),
  L("Sem critério: melhor RMSE do IC nas 16 comuns", 1.524, 0.0005, contraste$RMSE_22[contraste$variavel == "IC"]),
  L("Distância média à estação mais próxima (km)", 30.8, 0.05, resumo_dist$dist_media_km),
  L("Percentil 75 da distância (km)", 41.8, 0.05, resumo_dist$p75_km),
  L("Percentil 90 da distância (km)", 51.9, 0.05, resumo_dist$p90_km),
  L("Distância máxima (km)", 78.7, 0.05, resumo_dist$dist_max_km),
  L("Área a mais de 50 km de uma estação (%)", 12.3, 0.05, resumo_dist$pct_area_acima_50km),
  L("Municípios a mais de 40 km", 41, 0, leitura$municipios_acima_40km),
  L("IC na grade: mínimo (°C)", 29.2, 0.05, fx("IC", "minimo")),
  L("IC na grade: máximo (°C)", 36.0, 0.05, fx("IC", "maximo")),
  L("NET na grade: mínimo (°C)", 22.5, 0.05, fx("NET", "minimo")),
  L("NET na grade: máximo (°C)", 28.1, 0.05, fx("NET", "maximo")),
  L("Correlação IC x NET municipal", 0.918, 0.0005, leitura$r_IC_NET_municipal),
  L("Municípios com deslocamento > 30 posições", 14, 0, leitura$desloc_mais_30),
  L("Municípios com deslocamento > 20 posições", 52, 0, leitura$desloc_mais_20),
  L("r(deslocamento, V) municipal", -0.54, 0.005, leitura$r_deslocamento_V),
  L("Municípios NWS cuidado", 66, 0, cls("National Weather Service", "^Cuidado \\(", "municipios")),
  L("População NWS cuidado", 2086866, 0, cls("National Weather Service", "^Cuidado \\(", "populacao")),
  L("Municípios NWS cuidado extremo", 101, 0, cls("National Weather Service", "^Cuidado extremo", "municipios")),
  L("Municípios Fortaleza confortável", 72, 0, cls("Ajustada para Fortaleza", "^Confort", "municipios")),
  L("Municípios Fortaleza desconfortável", 95, 0, cls("Ajustada para Fortaleza", "^Desconfort", "municipios")),
  L("Municípios a menos de 1 RMSE de um limiar do NWS", 90, 0, leitura$perto_limiar_municipios),
  L("Municípios fora da zona de alta confiança", 38, 0, leitura$fora_zona_alta_municipios),
  L("População fora da zona de alta confiança", 275037, 0, leitura$fora_zona_alta_populacao),
  L("IC > 32 °C às 15h: menor frequência entre estações (%)", 6.5, 0.05, min(freq$pct_IC_acima_32)),
  L("IC > 32 °C às 15h: maior frequência entre estações (%)", 97.3, 0.05, max(freq$pct_IC_acima_32)),
  L("IC > 41 °C às 15h: maior frequência entre estações (%)", 0.9, 0.05, max(freq$pct_IC_acima_41)),
  L("IC > 41 °C às 15h: ocorrências", 5, 0, sum(freq$ocorrencias_acima_41)),
  L("Natal: dias comuns às 15h", 216, 0, if (!is.null(natal)) natal$dias_comuns[1] else NA),
  L("Natal: diferença do IC, INMET menos EMPARN (°C)", 1.55, 0.005,
    if (!is.null(natal)) natal$diferenca_INMET_menos_EMPARN[natal$variavel == "IC (°C)"] else NA)) %>%
  mutate(diferenca = obtido - publicado,
         situacao = dplyr::case_when(is.na(obtido) ~ "não calculado",
                                     abs(diferenca) <= tolerancia + 1e-9 ~ "confere",
                                     TRUE ~ "difere"))
gravar(conf, dO, "conferencia_valores_publicados.csv")
msg("   ", sum(conf$situacao == "confere"), " de ", nrow(conf), " valores conferem com o artigo")
print(as.data.frame(conf %>% filter(situacao != "confere") %>% select(item, publicado, obtido)), digits = 6)

# Comparação opcional com as predições originais da validação cruzada (por método)
for (h in P$horas) {
  arq_orig <- file.path(DIR_ORIGINAL, sprintf("loocv_pred_h%02d.csv", h))
  if (!existe(arq_orig)) next
  orig <- readr::read_csv(arq_orig, show_col_types = FALSE)
  if (!all(c("var", "metodo", "estacao", "pred") %in% names(orig))) next
  cmp <- cv %>% filter(horario == h) %>% left_join(METODOS %>% select(metodo = id, codigo), by = "metodo") %>%
    mutate(k = chave(estacao)) %>%
    inner_join(orig %>% transmute(variavel = var, codigo = metodo, k = chave(estacao), pred_orig = pred),
               by = c("variavel", "codigo", "k")) %>%
    group_by(variavel, metodo, codigo) %>%
    summarise(estacoes = n(), maior_diferenca = max(abs(pred - pred_orig), na.rm = TRUE), .groups = "drop")
  gravar(cmp, dO, sprintf("comparacao_predicoes_originais_h%02d.csv", h))
  msg(sprintf("   %02dh: %d de %d combinações de variável e método reproduzem as predições originais (diferença < 0,001)",
              h, sum(cmp$maior_diferenca < 1e-3, na.rm = TRUE), nrow(cmp)))
}
writeLines(capture.output(utils::sessionInfo()), file.path(dO, "sessao_R.txt"))


## ============================================================================
## Apêndice. Código do Google Earth Engine que gera as covariáveis
## ============================================================================
# Cole o arquivo gravado em 16_apendice/gee_covariaveis_RN_2025.js no Code Editor
# (https://code.earthengine.google.com) e execute as tarefas da aba Tasks. Os
# GeoTIFFs vão para a pasta covariaveis_RN do Google Drive; copie-os para DIR_COV.
dAp <- pasta("16_apendice")
writeLines(c(
"// Covariáveis do artigo: SRTM, NDVI (MOD13Q1 v061, média de 2025) e ERA5-Land",
"// (médias anuais às 03, 12 e 18 UTC dos índices calculados hora a hora).",
"var regiao = ee.Geometry.Rectangle([-38.9, -7.2, -34.7, -4.6], null, false);",
"var Z0 = 0.03;",
"var FATOR = Math.log(2 / Z0) / Math.log(10 / Z0);   // 0,723",
"",
"function indices(img) {",
"  var T  = img.select('temperature_2m').subtract(273.15);",
"  var Td = img.select('dewpoint_temperature_2m').subtract(273.15);",
"  // Umidade relativa pela fórmula de Magnus (Alduchov e Eskridge 1996)",
"  var es = T.multiply(17.625).divide(T.add(243.04)).exp();",
"  var ea = Td.multiply(17.625).divide(Td.add(243.04)).exp();",
"  var UR = ea.divide(es).multiply(100).min(100);",
"  var V10 = img.select('u_component_of_wind_10m').hypot(img.select('v_component_of_wind_10m'));",
"  var V2 = V10.multiply(FATOR);",
"  // IC: Rothfusz (1990) com os ajustes do NWS, em °F, convertido para °C",
"  var Tf = T.multiply(1.8).add(32);",
"  var simples = Tf.expression('0.5 * (Tf + 61 + (Tf - 68) * 1.2 + RH * 0.094)', {Tf: Tf, RH: UR});",
"  var cheio = Tf.expression(",
"    '-42.379 + 2.04901523 * Tf + 10.14333127 * RH - 0.22475541 * Tf * RH' +",
"    ' - 0.00683783 * Tf * Tf - 0.05481717 * RH * RH + 0.00122874 * Tf * Tf * RH' +",
"    ' + 0.00085282 * Tf * RH * RH - 0.00000199 * Tf * Tf * RH * RH', {Tf: Tf, RH: UR});",
"  var c1 = UR.lt(13).and(Tf.gte(80)).and(Tf.lte(112));",
"  var a1 = ee.Image(13).subtract(UR).divide(4).multiply(",
"    ee.Image(17).subtract(Tf.subtract(95).abs()).max(0).divide(17).sqrt()).multiply(-1);",
"  var c2 = UR.gt(85).and(Tf.gte(80)).and(Tf.lte(87));",
"  var a2 = UR.subtract(85).divide(10).multiply(ee.Image(87).subtract(Tf).divide(5));",
"  var ajuste = ee.Image(0).where(c1, a1).where(c2, a2);",
"  var hi = simples.where(simples.add(Tf).divide(2).gte(80), cheio.add(ajuste));",
"  var IC = hi.subtract(32).divide(1.8);",
"  // NET: Li e Chan (2000), vento a 2 m",
"  var den = ee.Image(0.68).subtract(UR.multiply(0.0014))",
"    .add(ee.Image(1).divide(V2.pow(0.75).multiply(1.4).add(1.76)));",
"  var NET = ee.Image(37).subtract(ee.Image(37).subtract(T).divide(den))",
"    .subtract(T.multiply(0.29).multiply(ee.Image(1).subtract(UR.divide(100))));",
"  return ee.Image.cat([T.rename('T_era5'), Td.rename('Td_era5'), UR.rename('UR_era5'),",
"    V10.rename('V10_era5'), V2.rename('V2m_era5'), IC.rename('IC_era5'), NET.rename('NET_era5_v2m')]);",
"}",
"",
"[3, 12, 18].forEach(function (h) {",
"  var media = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')",
"    .filterDate('2025-01-01', '2026-01-01')",
"    .filter(ee.Filter.calendarRange(h, h, 'hour'))",
"    .map(indices).mean();",
"  var nome = 'RN_era5land_h' + (h < 10 ? '0' + h : h) + '_2025';",
"  Export.image.toDrive({image: media.toFloat(), description: nome, fileNamePrefix: nome,",
"    folder: 'covariaveis_RN', region: regiao, scale: 11132, crs: 'EPSG:4326', maxPixels: 1e10});",
"});",
"",
"var ndvi = ee.ImageCollection('MODIS/061/MOD13Q1').filterDate('2025-01-01', '2026-01-01')",
"  .select('NDVI').mean().multiply(0.0001).rename('ndvi');",
"Export.image.toDrive({image: ndvi.toFloat(), description: 'RN_ndvi_2025', fileNamePrefix: 'RN_ndvi_2025',",
"  folder: 'covariaveis_RN', region: regiao, scale: 250, crs: 'EPSG:4326', maxPixels: 1e10});",
"",
"var srtm = ee.Image('USGS/SRTMGL1_003').select('elevation').rename('altitude');",
"Export.image.toDrive({image: srtm.toFloat(), description: 'RN_srtm', fileNamePrefix: 'RN_srtm',",
"  folder: 'covariaveis_RN', region: regiao, scale: 90, crs: 'EPSG:4326', maxPixels: 1e10});",
"",
"// A LCZ (Demuzere et al. 2022) foi obtida com o pacote LCZ4r no R",
"// (BAIXAR_LCZ_COM_LCZ4R = TRUE no script)."),
file.path(dAp, "gee_covariaveis_RN_2025.js"))

msg("Concluído. Resultados em ", DIR_SAIDA)
