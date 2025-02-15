# **spnhppzi: Modelagem Bayesiana de Eventos Recorrentes com Inflação de Zeros e Correlação Espacial**  

O **spnhppzi** é um pacote desenvolvido para a **modelagem de dados de eventos recorrentes** que apresentam **inflação de zeros e correlação espacial**. Ele adota uma **abordagem bayesiana**, implementando **modelos hierárquicos** para capturar **estruturas complexas** associadas à repetição de eventos e suas dependências espaciais.  

A modelagem espacial é baseada no **modelo Condicional Intrinsecamente Autorregressivo (ICAR)**, que permite incorporar a correlação espacial de maneira eficiente. Além disso, o pacote oferece **modelos paramétricos e semiparamétricos**, nos quais **polinômios de Bernstein** são empregados para a **modelagem da função de intensidade de linha de base**. Essa abordagem aumenta a flexibilidade do modelo, tornando-o aplicável a cenários onde **modelos estritamente paramétricos podem falhar na representação da complexidade dos dados**.  

## **Simulação de Dados**  

O **spnhppzi** também inclui **funções para simulação de dados de eventos recorrentes**, utilizando uma adaptação do pacote **SIMREC** ([Farrington et al., 2014](https://doi.org/10.18637/jss.v058.i02)). Isso permite a geração de **dados espacialmente correlacionados**, possibilitando a **avaliação do desempenho dos modelos** e experimentação com diferentes cenários de recorrência e dependência espacial.  

## **Referências e Aplicação**  

Este repositório contém as funções utilizadas para gerar os resultados apresentados no artigo:  

📄 **"The Analysis of Criminal Recidivism: A Hierarchical Model-Based Approach for the Analysis of Zero-Inflated, Spatially Correlated Recurrent Events Data"**, disponível em [arXiv:2405.02666](https://arxiv.org/abs/2405.02666).  

Detalhes adicionais sobre a metodologia podem ser encontrados na tese:  

📖 **"Modelos Hierárquicos para a Análise de Dados de Eventos Recorrentes com Inflação de Zeros e Correlação Espacial"**, disponível mediante solicitação ao autor pelo e-mail **alisson.ccs2@gmail.com**.  

## **Funcionalidades**  

- **Modelagem hierárquica bayesiana** para dados de eventos recorrentes;  
- **Modelagem da correlação espacial** entre eventos por meio do **modelo ICAR**;  
- **Versões paramétricas e semiparamétricas**, com **polinômios de Bernstein** para maior flexibilidade na modelagem da função de intensidade de linha de base;  
- **Simulação de dados** de eventos recorrentes espacialmente correlacionados (adaptação do **SIMREC**);  
- **Ferramentas para pré-processamento de dados e estimação de modelos**.  

## **Contribuição e Contato**  

Agradecemos seu interesse no **spnhppzi**! Caso encontre algum problema ou tenha sugestões, fique à vontade para **abrir uma issue** ou entrar em contato pelo e-mail **alisson.ccs2@gmail.com**.  
