# Synthetic EMR Data

This synthetic dataset was created for verifying the pyPheWAS package. We have made it freely available to allow users to test pyPheWAS's capabilities.

### Demographic Summary

The dataset contains 10,000 subjects split evenly between cases (Dx=1)
and controls (Dx=0), and includes two demographic variables: biological
sex and maximum age at visit.

|                | Subjects | Sex [% Female] | Max Age At Visit [mean (std.)] |
| -------------- |:-------: | :------------: | :----------------------------: |
| Case (Dx=1)    | 5,000    | 70%            | 59.946 (9.563)                 |
| Control (Dx=0) | 5,000    | 40%            | 60.802 (9.448)                 |

### ICD Record Generation
A mix of 103,493 ICD-9 and ICD-10 code events were generated for this dataset,
covering 31 PheCode associations. Three types of PheCode associations were created, including:
* 20 **background** associations,
* 9 **primary** associations,
* and 2 **confounding** associations.

**Background** | insignificant associations between Dx and the PheCode

ICD events were generated such that each background PheCode would have a small pre-specified effect size, randomly generated via a uniform distribution over the range [-0.1, 0.1]; individuals’ ages for each event were randomly generated using a uniform distribution over the range [30, 50]. *pyPheWAS should accurately estimate each background association’s effect size but determine that the association is insignificant.*

**Primary** | true associations between Dx and the PheCode

ICD events were generated such that primary PheCodes would have a unique pre-specified effect size (log odds ratio) across the full cohort; individuals’ ages for each event were randomly generated using a uniform distribution over the range [30, 50]. The nine primary PheCodes and their respective effect
sizes are shown in the table below.

| PheCode | Phenotype                                        | Log Odds Ratio |
| ------- |------------------------------------------------- | -------------- |
| 338.2   | Chronic pain                                     | 1.50           |
| 340     | Migraine                                         | 1.10           |
| 1011    | Complications of surgical and medical procedures | 0.70           |
| 296.22  | Major depressive disorder                        | 0.60           |
| 530.11  | GERD                                             | 0.30           |
| 401     | Hypertension                                     | 0.25           |
| 041     | Bacterial infection NOS                          | -0.20          |
| 1009    | Injury, NOS                                      | -0.60          |
| 495     | Asthma                                           | -1.00          |

*pyPheWAS should accurately estimate each primary association’s effect size and determine that the association is statistically significant.*


**Confounding** | false positive associations caused by the confounding effect of either sex or age

PheCode 174.1 (Breast cancer [female]) was used as a **sex-confounded** PheCode. ICD events were generated such that all females in the dataset had equal odds of having PheCode 174.1 in their record; ages for each event were randomly generated using a uniform distribution over the range [30, 50]. Because females were disproportionally represented across the case and control groups, however, the PheCode’s cohort-wide effect size is positively skewed to a 0.6 log odds ratio.

PheCode 292.2 (Mild cognitive impairment) was used as an **age-confounded** PheCode. ICD events were generated such that PheCode 292.2 would have a -0.2 log odds ratio; however, event ages were randomly generated using a uniform distribution over the higher age range [65,70]. This resulted in PheCode 292.2 being highly associated with larger values of MAV.

*Without controlling for the confounding variable, pyPheWAS should identify a significant association with these confounded PheCodes; including the confounding variable as a covariate, however, should reduce (or eliminate) the confounded association.*
