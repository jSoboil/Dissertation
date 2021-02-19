# A Sensitivity Analysis Framework for Health Economic Evaluation in Middle Income Countries: Appropriately Incorporating a Comprehensive Approach

<img src="misc/logo.jpg" width="260" align="right" />
<br/>

![CEAbadge](https://img.shields.io/github/issues/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/last-commit/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/license/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/badge/R-v4.0.0+-blue)
![CEAbadge](https://img.shields.io/badge/JAGS-v4.3.0-blue)

## Authors
1. [Joshua Soboil](https://orcid.org/0000-0003-1362-8465)<sup>a,b</sup>
2. [Lucy Cunnama](https://orcid.org/0000-0003-2134-4905)<sup>b</sup>
3. [Tommy Wilkinson](https://orcid.org/0000-0003-0806-2196)<sup>c</sup>

<sup>a. Corresponding author's [email](mailto:soboil.joshua@gmail.com) <br/>
b. Health Economics Unit, School of Public Health and Family Medicine, University of Cape Town <br/>
c. World Bank Group, 1818 H Street Washington, DC <br/>
<sup>
<br/>

## Brief
This repository stores the code to run a hypothetical case study assessing the practical value of a comprehensive approach to health economic decision-modelling in Middle-Income Country contexts. The model is coded in the R language and is a recreated model originally developed by [Sinanovic E. et al.](https://doi.org/10.1016/j.vaccine.2009.08.004), titled:

>The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical cancer screening programme in South Africa.

<br/>

## Abstract
Estimating uncertainty provides indispensable information on a decision’s set of possible outcomes. Within health economics specifically, there are several conventional sensitivity analysis methods to choose from, each with unique advantages and disadvantages to consider. Accordingly, when constructing a health economic decision model, it is critical to select a sensitivity analysis approach appropriate for the given decision context. This point is particularly salient to Middle-Income Countries, where there is relatively heightened resource scarcity and thus increased opportunity-costs. Many Middle-Income Countries face an acute shortage of accessible as well as reliable evidence, resulting in a frequent practice of imputing data derived from external jurisdictions. Conversely, there is frequent shortages of technical skills as well as research capacity and thus a strong complementary need to consider the contextual feasibility of applying resource demanding methodologies. Given the availability of several sensitivity analysis methods which can sufficiently accommodate diverse areas of uncertainty, it is critical to establish whether and when technical benefits of more complex methods results in real-world value. In light of the above, critical insight is gained through a hypothetical case study which applies a comprehensive approach to decision-modelling, implemented using the R and JAGS languages. The case study replicates a decision model originally used to inform the potential cost-effectiveness of adding a bivalent Human Papilloma Virus vaccine to South Africa’s primary health care cervical cancer screening programme. Importantly however, the case study provides a practical comparison of potential pros and cons of implementing more complex techniques and, perhaps most significantly, guides the establishment of a sensitivity analysis selection framework guided by a conceptual ‘fruitfulness’ (i.e. a bang-for-buck) centred on several fundamental decision-making categories. By establishing this framework, we argue that it encourages Middle-Income Country analysts to select sensitivity analysis methods more judiciously, will help reduce the methodological variation apparent in Middle-Income Country settings, and simultaneously provide decision-makers with greater methodological transparency in the selection of sensitivity analysis methods.

<details>
<summary>Technical note</summary>
Before running the model, ensure that the local working directory is set to the location of the .Rproj folder saved on your computer. In RStudio, the easiest way to select the local directory path is by pressing Ctrl + Shift + H.

The coding style throughout the model follows the framework proposed by [Alarid-Escudero F. et al.](https://doi.org/10.1007/s40273-019-00837-x) titled:

>A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling. 

<br/>