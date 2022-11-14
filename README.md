# A Sensitivity Analysis Framework for Health Economic Evaluation in Middle Income Countries: Appropriately Incorporating a Comprehensive Approach

<img src="misc/logo.jpg" width="260" align="right" />
<br/>

![CEAbadge](https://img.shields.io/github/issues/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/last-commit/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/github/license/jSoboil/Dissertation)
![CEAbadge](https://img.shields.io/badge/R-v4.0.0+-blue)
![CEAbadge](https://img.shields.io/badge/JAGS-v4.3.0-blue)

## Authors
1. [Joshua Soboil, MPH](https://orcid.org/0000-0003-1362-8465)<sup>a,b</sup>
2. [Lucy Cunnama, PhD](https://orcid.org/0000-0003-2134-4905)<sup>b</sup>
3. [Tommy Wilkinson, MSc](https://orcid.org/0000-0003-0806-2196)<sup>c</sup>

<sup>a. Corresponding author's [email](mailto:joshua@soboils.com?subject=[GitHub]%20Dissertaion%20Enquiry) <br/>
b. Health Economics Unit, School of Public Health and Family Medicine, University of Cape Town <br/>
c. World Bank Group <br/>
<sup>
<br/>

## Technical Brief
This repository stores the code to run a hypothetical model which is used to assess the practical value of a comprehensive approach to health economic decision-modelling in Middle-Income Country contexts. The replicated model is primarily coded in the R language, but it is supplemented with JAGS. Note that the replicated model is based on an original model developed by [Sinanovic E. et al.](https://doi.org/10.1016/j.vaccine.2009.08.004), titled:

>The potential cost-effectiveness of adding a human papillomavirus vaccine to the cervical cancer screening programme in South Africa`

<br/>

Before running the replicated model, please ensure that the local working directory is set to the location of the .Rproj folder on your computer. In RStudio, the easiest way to select the local directory path is by pressing Ctrl + Shift + H. You can also download JAGS [here](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/) if necessary. The coding style throughout follows the framework proposed by [Alarid-Escudero F. et al.](https://doi.org/10.1007/s40273-019-00837-x) titled:

>A Need for Change! A Coding Framework for Improving Transparency in Decision Modeling. 

<br/>

## Abstract
When constructing a health economic decision model it is critical to select a sensitivity analysis approach that is appropriate for the decision context. This point is particularly salient to Middle-Income Countries (MICs), where there is relatively heightened resource scarcity and increased opportunity-cost. MICs face acute shortages of accessible as well as high-quality evidence, resulting in frequent imputing of data from external jurisdictions. Conversely, there are shortages in skills and research capacity, creating a strong complementary need to consider the contextual feasibility of applying resource demanding methodologies. Given the above, it is critical to establish when the technical benefits of more complex and resource demanding sensitivity analysis methods result in real-world value. We apply a comparative case study using a comprehensive approach to decision-modelling, implemented in the R and JAGS languages. The case study replicates a deterministic model originally used to inform the cost-effectiveness of adding a bivalent Human Papilloma Virus (HPV) vaccine to South Africa’s public health care cervical cancer screening programme. Crucially, the case study provides critical insight into the potential pros and cons of implementing more complex sensitivity analysis techniques within MIC climates. Our findings indicate that the benefits of more advanced sensitivity analysis methods are nuanced; are contextually beneficial according to a case-by-case basis; and, moreover, that choosing a sensitivity analysis method should be guided by a conceptual ‘fruitfulness’ (i.e., a bang-for-buck). To aid analysts in the process of selecting an appropriate approach to achieving a sensitivity analysis output using the insight gained from the comparative case study, we provide a general framework containing three core conceptual areas, namely: Decision-Maker Preferences (Investment, Decision Power, and Risk Aversion), Analytical Considerations (Available resources and Indirect Evidence) and Policy Context (Knowledge of Topic and Technical Expertise). The framework intends to encourage a judicious selection of sensitivity analysis methods, reduce the methodological variation apparent in MIC settings, and provide health care decision-makers with greater methodological transparency.