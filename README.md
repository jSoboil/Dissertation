---
title: A Selection Framework for Incorporating a Fully Integrated Bayesian Approach to Cost-Effectiveness Evaluation in Middle Income Country Contexts
css: misc/css/styles.css
---

<br/>

<img src="misc/logo.jpg" width="260" align="right" />

<br/>

![CEAbadge](https://img.shields.io/github/issues/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/last-commit/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/license/jSoboil/Dissertation?label=license)
![CEAbadge](https://img.shields.io/badge/R-v4.0.0+-blue)
![CEAbadge](https://img.shields.io/badge/JAGS-v4.3.0-blue)

<br/>

## Authors
1. [Joshua Soboil](https://www.linkedin.com/in/joshua-soboil-067351172/)<sup>a,b</sup>
2. [Lucy Cunnama](https://scholar.google.co.za/citations?hl=en&user=eG7OJ7EAAAAJ)<sup>b</sup>
3. [Tommy Wilkinson](https://twitter.com/Tommy_HealthSA)<sup>b</sup>

<sup>a. Corresponding author ([email](mailto:soboil.joshua@gmail.com)) \
b. Health Economics Unit, School of Public Health and Family Medicine, University of Cape Town.<sup>

<br/>

## Brief
<p>This repository stores a Cost-Effectiveness Analysis model coded in the R language. The model is is a replication of an original model developed by [Sinanovic E. et al.](https://doi.org/10.1016/j.vaccine.2009.08.004), titled:

>The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical cancer screening programme in South Africa.


## Background
<p>Due to the scarcity of resources faced by MICs, it is imperative that economic evaluation methods within such contexts are selected judiciously. Accordingly, the article associated with this repository proposes a selection framework intended to streamline and improve the choice of sensitivity analysis methods for health economic evaluations conducted in MICs.<p>

<p>The framework was established by comparing a MIC cost-effectiveness analysis model of the human papillomavirus (HPV) bivalent vaccine to a recreated version using more comprehensive methods and software. The primary intention of the replication was to establish whether a more comprehensive approach would automatically result in higher quality decision-making. The original study was developed in TreeAge, employed a deterministic approach to sensitivity analysis, and found the intervention to be more costly and more effective (while citing considerable uncertainty surrounding the cost of the vaccine). In contrast, the replication study chose a more integrated approach to decision modelling, using Bayesian ‘Markov Chain Monte Carlo’ (MCMC). The replication exercise indicated that the original model under-estimated decision uncertainty however, given the relatively low costs of the intervention, this was negligible and did not result in any significant change to the overall cost-effectiveness of the intervention. Thus, the authors argue that the original model was sufficient in its quality as well as complexity to answer the decision-problem.<p>


## Discussion
<p>The result of the replication study suggests that the *technical* benefits of more complex methods do not always translate into *practical* benefits. Moreover, because social values play a fundamental role in informing health care resource allocation decisions, it is imperative that the analyst always consider the perspective of the decision-maker. Health economic decision-models are not solely predictive tools, but also communication tools that help decision-makers to weigh both the social and economic consequences of their resource allocation decisions. As such, the authors propose a selection framework using the Pugh-matrix concept selection method. Several categories associated with decision-maker preferences, available resources, as well as the type of evidence and software available are defined as essential to consider when choosing the appropriate sensitivity analysis method within MIC contexts. However, because of its concrete and qualitative nature, the authors anticipate that the selection framework will additionally serve as a beneficial instrument for better knowledge translation between an analyst and decision-maker. See a rough version of the selection matrix below<p>
<br/>

<img src="figs/Pugh_matrix.png" width="490" style="float: centre" />

<br/>
<br/>

<details>
<summary>Technical note</summary>
<p>Before running the model, ensure that the local working directory is set to the location of the .Rproj folder saved on your computer. In RStudio, the easiest way to select the local directory path is by pressing Ctrl + Shift + H.<p>

<p>The coding style throughout the model follows the framework proposed by [Alarid-Escudero F. et al.](https://doi.org/10.1007/s40273-019-00837-x) titled:

>A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling. 

<br/>