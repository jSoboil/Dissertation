---
title: A Selection Framework for Incorporating a Fully Integrated Bayesian Approach to Cost-Effectiveness Evaluation in Middle Income Country Contexts
output:
  md_document: 
css: "misc/css/styles.css"
---
\
\
\
![](misc/logo.jpg){width=25% align=right}
\
\
\
\
\
\
![CEAbuild](https://img.shields.io/github/issues/jSoboil/Dissertation)
\
\
\
\

## Background
<p>This repository stores the replicated Cost-Effectiveness Analysis code used to explore the strengths of a fully integrated Bayesian Posterior simulation approach to Cost-Effectiveness Analysis in a middle-income country context. The [original study](https://www.sciencedirect.com/science/article/pii/S0264410X09011670?via%3Dihub) used Deterministic Sensitivity Analysis (DSA) and was developed by Sinanovic E., et al. (2009) titled: 

>The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical cancer screening programme in South Africa.

<p>Due to the scarcity of resources faced by MICs, it is imperative that economic evaluation methods within such contexts are selected judiciously. Accordingly, the following paper proposes a selection framework intended to streamline and improve the choice of sensitivity analysis methods for health economic evaluations conducted in MICs.<p>

<p>The framework was established by comparing a MIC cost-effectiveness analysis model of the human papillomavirus (HPV) bivalent vaccine to a recreated version using more comprehensive methods and software. The primary intention of the replication was to establish whether a more comprehensive approach would automatically result in higher quality decision-making. The original study was developed in TreeAge, employed a deterministic approach to sensitivity analysis, and found the intervention to be more costly and more effective (while citing considerable uncertainty surrounding the cost of the vaccine). In contrast, the replication study chose a more integrated approach to decision modelling, using Bayesian ‘Markov Chain Monte Carlo’ (MCMC). The replication exercise indicated that the original model under-estimated decision uncertainty however, given the relatively low costs of the intervention, this was negligible and did not result in any significant change to the overall cost-effectiveness of the intervention. Thus, the authors argue that the original model was sufficient in its quality as well as complexity to answer the decision-problem.<p>

<p>Accordingly, the findings of the replication study suggest that the technical benefits of more complex methods do not always translate into practical benefits. Moreover, because social values play a fundamental role in informing health care resource allocation decisions, it is imperative that the analyst always consider the perspective of the decision-maker. Health economic decision-models are not solely predictive tools, but also communication tools that help decision-makers to weigh both the social and economic consequences of their resource allocation decisions. As such, the authors propose a selection framework using the Pugh-matrix concept selection method. Several categories associated with decision-maker preferences, available resources, as well as the type of evidence and software available are defined as essential to consider when choosing the appropriate sensitivity analysis method within MIC contexts. However, because of its concrete and qualitative nature, the authors anticipate that the selection framework will additionally serve as a beneficial instrument for better knowledge translation between an analyst and decision-maker. See a rough version of the selection matrix below<p>
\
\
\
\
\
![](figs/Pugh_matrix.png){width=25% height=450px align=center style="float"}
\
\
\
\
\
\
\

## Technical Note
<p>Before running the model, ensure that the working directory is set to the location of the .Rproj folder when in a local directory. In RStudio, the easiest way to do this is by pressing Ctrl + Shift + H to select the local directory path.<p>
<details>
<summary>Coding style</summary>
<p>The coding styled used throughout this cost-effectiveness model follows the coding framework proposed by Alarid-Escudero et al. (2019) titled:

>A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling. 

Access the article [here](https://doi.org/10.1007/s40273-019-00837-x)<p>

