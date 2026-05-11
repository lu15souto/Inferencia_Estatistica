###############################################################
# Instalando e carregando pacotes
###############################################################

pacotes <- c(
  "quantmod", "PerformanceAnalytics", "dplyr", "tidyr", "ggplot2", "ggpubr",
  "qqplotr", "car", "lawstat", "multcompView", "PMCMRplus", "reshape2",
  "moments", "tseries"
)

for(p in pacotes){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}


###############################################################
# Baixar preço das ações
###############################################################

tickers <- c("VALE3.SA","ITUB4.SA","PETR4.SA","BBDC4.SA","SBSP3.SA","BPAC11.SA",
             "ITSA4.SA","WEGE3.SA","ABEV3.SA","RDOR3.SA","EQTL3.SA","VIVT3.SA")

inicio <- "2000-01-01"
fim     <- "2025-11-26"

precos <- list()
for (t in tickers){
  try({
    dados <- getSymbols(t, src="yahoo", from=inicio, to=fim, auto.assign=FALSE)
    precos[[t]] <- Ad(dados)
  })
}

precos_xts <- do.call(merge, precos)
retornos <- na.omit(ROC(precos_xts, type="discrete"))


###############################################################
# Transformar em data frame longo
###############################################################

ret_df <- data.frame(Data = index(retornos),
                     coredata(retornos)) %>%
  pivot_longer(cols = -Data,
               names_to="Ativo",
               values_to="Retorno") %>%
  drop_na()

# garantir fator limpo
ret_df$Ativo <- droplevels(factor(ret_df$Ativo))


###############################################################
# Testes
###############################################################

# Normalidade (Shapiro–Wilk)
shapiro_results <- ret_df %>%
  group_by(Ativo) %>%
  summarise(
    p_value = shapiro.test(Retorno)$p.value,
    Normalidade = ifelse(p_value > 0.05, "Normal", "Não-normal")
  )
print(shapiro_results)

# Homogeneidade (Levene)
levene <- levene.test(ret_df$Retorno, ret_df$Ativo)
print(levene)

###############################################################
# Kruskal-Wallis
###############################################################

kw <- kruskal.test(Retorno ~ Ativo, data = ret_df)
print(kw)


###############################################################
#  post-hoc Dunn — sem ajuste e com ajuste de Bonferroni
###############################################################

ret_df <- ret_df %>% filter(is.finite(Retorno))
ret_df$Ativo <- droplevels(factor(ret_df$Ativo))


# Completendo a matriz, qe sera gerada no post-hoc Dunn, deixando a diagonal principal com 1, para facilciar a plotagem e a nalise dos resultados
complete_pmatrix <- function(pmat){
  ativos <- union(rownames(pmat), colnames(pmat))
  full <- matrix(1, nrow = length(ativos), ncol = length(ativos))
  rownames(full) <- ativos
  colnames(full) <- ativos
  
  for(i in 1:nrow(pmat)){
    for(j in 1:ncol(pmat)){
      full[rownames(pmat)[i], colnames(pmat)[j]] <- pmat[i,j]
      full[colnames(pmat)[j], rownames(pmat)[i]] <- pmat[i,j]
    }
  }
  return(full)
}

# --- Dunn sem ajuste
dunn_none <- PMCMRplus::kwAllPairsDunnTest(
  Retorno ~ Ativo,
  data = ret_df,
  p.adjust.method = "none"
)
pmat_none <- dunn_none$p.value
pmat_none_full <- complete_pmatrix(pmat_none)

# Definindo a diagonal principal 1
diag(pmat_none_full) <- 1

# Gerando letras de agrupamento
letras_none <- multcompView::multcompLetters(pmat_none_full)$Letters

# Dunn Bonferron
dunn_bonf <- PMCMRplus::kwAllPairsDunnTest(
  Retorno ~ Ativo,
  data = ret_df,
  p.adjust.method = "bonferroni"
)
pmat_bonf <- dunn_bonf$p.value
pmat_bonf_full <- complete_pmatrix(pmat_bonf)

# Definindo a diagonal principal 1
diag(pmat_bonf_full) <- 1

# Gerando letras de agrupamento
letras_bonf <- multcompView::multcompLetters(pmat_bonf_full)$Letters


###############################################################
# Criando tabelas dos agrupamentos
###############################################################

tabela_none <- data.frame(
  Ajuste = "Sem Ajuste",
  Ativo = names(letras_none),
  Grupo = letras_none
)

tabela_bonf <- data.frame(
  Ajuste = "Bonferroni",
  Ativo = names(letras_bonf),
  Grupo = letras_bonf
)

print(tabela_none)
print(tabela_bonf)



###############################################################
# 8) Gráfico dos comparativos dos grupos de letars
###############################################################

medias <- ret_df %>%
  group_by(Ativo) %>%
  summarise(
    media = mean(Retorno),
    sd = sd(Retorno),
    n = n(),
    se = sd/sqrt(n),
    ci = qt(0.975, df=n-1)*se
  )

# Sem ajuste
df_none_plot <- medias %>% left_join(tabela_none, by="Ativo")

ggplot(df_none_plot, aes(x=Ativo, y=media)) +
  geom_point(size=3, color="red") +
  geom_errorbar(aes(ymin=media-ci, ymax=media+ci), width=0.15) +
  geom_text(aes(label=Grupo, y=media+ci*1.1), size=5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Dunn Sem Ajuste", y="Retorno Médio")

# com ajuste de Bonferron
df_bonf_plot <- medias %>% left_join(tabela_bonf, by="Ativo")

ggplot(df_bonf_plot, aes(x=Ativo, y=media)) +
  geom_point(size=3, color="blue") +
  geom_errorbar(aes(ymin=media-ci, ymax=media+ci), width=0.15) +
  geom_text(aes(label=Grupo, y=media+ci*1.1), size=5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Dunn Bonferroni", y="Retorno Médio")


###############################################################
# Mapa de calor do p_valor
###############################################################

# Sem ajuste
dfh_none <- reshape2::melt(pmat_none_full)
colnames(dfh_none) <- c("Ativo1","Ativo2","p")

ggplot(dfh_none, aes(Ativo1, Ativo2, fill=p)) +
  geom_tile() +
  scale_fill_gradient(low="red", high="white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Heatmap p-values — Sem Ajuste")

#  Bonferroni
dfh_bonf <- reshape2::melt(pmat_bonf_full)
colnames(dfh_bonf) <- c("Ativo1","Ativo2","p")

ggplot(dfh_bonf, aes(Ativo1, Ativo2, fill=p)) +
  geom_tile() +
  scale_fill_gradient(low="red", high="white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Heatmap p-values — Bonferroni")

###############################################################
# Retorno acumulado
###############################################################

ret_acum <- cumprod(1 + retornos) - 1

df_acum <- data.frame(Data = index(ret_acum), coredata(ret_acum)) %>%
  pivot_longer(cols=-Data, names_to="Ativo", values_to="RetAcum")

ggplot(df_acum, aes(Data, RetAcum, color=Ativo)) +
  geom_line() +
  theme_bw() +
  labs(title="Retorno Acumulado (2020–2025)", y="Retorno Acumulado")

###############################################################
# Bloxplot da districuição dos retornos diários
###############################################################

ggplot(ret_df, aes(x=Ativo, y=Retorno, fill=Ativo)) +
  geom_boxplot(alpha=0.7) +
  geom_jitter(width=0.1, alpha=0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none") +
  labs(title="Distribuição dos Retornos Diários")

###############################################################
# Assimetria e Curtose
###############################################################

momentos <- ret_df %>%
  group_by(Ativo) %>%
  summarise(
    Assimetria = skewness(Retorno),
    Curtose = kurtosis(Retorno)
  )

print(momentos)

###############################################################
# Histograma da ditricuição de retorno
###############################################################

ggplot(ret_df, aes(Retorno, fill=Ativo)) +
  geom_histogram(aes(y=..density..), bins=60, alpha=0.6) +
  geom_density(color="black") +
  facet_wrap(~Ativo, scales="free") +
  theme_bw() +
  labs(title="Distribuição dos Retornos — Histograma + Densidade")

###############################################################
# Mapa de calor da correlação entre os ativos
###############################################################

cor_mat <- cor(retornos, use="pairwise")

df_cor <- reshape2::melt(cor_mat)
colnames(df_cor) <- c("Ativo1","Ativo2","Cor")

ggplot(df_cor, aes(Ativo1, Ativo2, fill=Cor)) +
  geom_tile() +
  scale_fill_gradient2(low="red", high="blue", mid="white", midpoint=0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Correlação dos Retornos")

###############################################################
# Risco X Retorno
###############################################################

risk_return <- retornos %>%
  apply(2, function(x) c(mean=mean(x), sd=sd(x))) %>%
  t() %>%
  as.data.frame()

risk_return$Ativo <- rownames(risk_return)

ggplot(risk_return, aes(x=sd, y=mean, label=Ativo)) +
  geom_point(size=3, color="blue") +
  geom_text(nudge_y=0.0005) +
  theme_bw() +
  labs(title="Retorno Médio vs Volatilidade",
       x="Desvio-Padrão (Risco)", y="Retorno Médio")

