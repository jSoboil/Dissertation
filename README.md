# A Sensitivity Analysis Selection Framework for Health Economic Evaluation in Middle Income Countries: Appropriately Incorporating a Comprehensive Approach

<img src="misc/logo.jpg" width="260" align="right" />
<br/>

![CEAbadge](https://img.shields.io/github/issues/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/last-commit/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/license/jSoboil/Dissertation?label=license)
![CEAbadge](https://img.shields.io/badge/R-v4.0.0+-blue)
![CEAbadge](https://img.shields.io/badge/JAGS-v4.3.0-blue)

## Authors
1. [Joshua Soboil](https://orcid.org/0000-0003-1362-8465)<sup>a,b</sup>
2. [Lucy Cunnama](https://orcid.org/0000-0003-2134-4905)<sup>b</sup>
3. [Tommy Wilkinson](https://orcid.org/0000-0003-0806-2196)<sup>b</sup>

<sup>a. Corresponding author's [email](mailto:soboil.joshua@gmail.com) <br/>
b. Health Economics Unit, School of Public Health and Family Medicine, University of Cape Town<sup>
<br/>

## Abstract
Sensitivity Analysis forms a vital part of health economic evaluation. Estimating the uncertainty of a decision can provide indispensable information on a decision’s set of potential outcomes. Moreover, there is a diverse set of established and developing approaches to choose from when conducting a sensitivity analysis, each associated with unique rewards and shortfalls that are important to consider. However, as is noted across health economic literature, it is imperative to select these methods judiciously, primarily by taking into account the context of a given decision problem. This is particularly pertinent to the context of Middle-Income Countries, as there is heightened scarcity of resources and thus increased opportunity-costs relative to higher-income jurisdictions. Likewise, many Middle-Income Countries face an acute shortage of accessible and reliable evidence which has resulted in the common practice of imputing data, into health economic decision models, from external jurisdictions. Given the general context of Middle-Income Countries, it is thus crucial to establish whether the technical benefits of more complex methods inevitably translate into tangible, practical benefits. Accordingly, by gaining critical insight from a hypothetical case study, we establish a framework guided by a conceptual ‘fruitfulness’ (i.e. a bang-for-buck), centred on several fundamental decision-making categories. We argue that the framework will encourage Middle-Income Country analysts to select sensitivity analysis methods more judiciously, help reduce the methodological variation apparent in Middle-Income Country settings, and simultaneously provide decision-makers with greater methodological transparency.

## Brief
This repository stores the code to run a hypothetical case study assessing the practical value of a comprehensive approach to health economic decision-modelling in Middle-Income Country contexts. The model is coded in the R language and is a recreated model originally developed by [Sinanovic E. et al.](https://doi.org/10.1016/j.vaccine.2009.08.004), titled:

>The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical cancer screening programme in South Africa.

<details>
<summary>Technical note</summary>
Before running the model, ensure that the local working directory is set to the location of the .Rproj folder saved on your computer. In RStudio, the easiest way to select the local directory path is by pressing Ctrl + Shift + H.

The coding style throughout the model follows the framework proposed by [Alarid-Escudero F. et al.](https://doi.org/10.1007/s40273-019-00837-x) titled:

>A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling. 

<br/>