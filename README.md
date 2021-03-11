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
c. World Bank Group <br/>
<sup>
<br/>

## Technical Brief
This repository stores the code to run a hypothetical case study assessing the practical value of a comprehensive approach to health economic decision-modelling in Middle-Income Country contexts. The model is coded in the R language and is based on an original model developed by [Sinanovic E. et al.](https://doi.org/10.1016/j.vaccine.2009.08.004), titled:

>The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical cancer screening programme in South Africa`

<br/>

Before running the model, ensure that the local working directory is set to the location of the .Rproj folder saved on your computer. In RStudio, the easiest way to select the local directory path is by pressing Ctrl + Shift + H. Please download JAGS [here](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/). The coding style throughout the model follows the framework proposed by [Alarid-Escudero F. et al.](https://doi.org/10.1007/s40273-019-00837-x) titled:

>A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling. 

<br/>

## Abstract
When constructing a health economic decision model, it is critical to select a sensitivity analysis approach appropriate for the decision context. This point is particularly salient to Middle-Income Countries (MICs), where there is relatively heightened resource scarcity and increased opportunity-cost. MICs face acute shortages of accessible as well as high-quality evidence, resulting in a frequent imputing of data from external jurisdictions. Conversely, there are also shortages in skills and research capacity, creating a strong complementary need to consider the contextual feasibility of applying more resource demanding sensitivity analysis methodologies.. Given the above, it is therefore critical to establish whether and when the technical benefits of complex and resource demanding methods result in real-world value. We apply a comparative case study using a comprehensive approach to decision-modelling, implemented in the R language and JAGS package. Specifically, the case study replicates a deterministic model originally used to inform the cost-effectiveness of adding a bivalent Human Papilloma Virus (HPV) vaccine to South Africa’s public health care cervical cancer screening programme. Crucially, the case study provides critical insight into the pros and cons of implementing more complex sensitivity analysis techniques within MIC climates. Our findings indicate that the benefits of more advanced sensitivity analysis methods are nuanced; are therefore contextually beneficial according to a case-by-case basis; and, moreover, choosing a sensitivity analysis method should be guided by a conceptual ‘fruitfulness’ (i.e. a bang-for-buck), more than a mere desire to reduce model complexity. To aid analysts in this process, from our comparative case study we provide a framework with three core concept areas namely Decision-Maker Preferences (Decision Power, Investment, Risk Aversion), Analytical Considerations (Available resources, Indirect Evidence) and Policy Context. The framework centred on these fundamental decision-making categories intends to encourage more judicious selection of sensitivity analysis methods; help reduce the methodological variation apparent in MIC settings; and simultaneously provide decision-makers with greater methodological transparency in the selection of sensitivity analysis methods. 
