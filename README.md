# Soil-Carbon-Debt

R code and spatial predictions for WHRC-TNC project modeling spatial extent of soil carbon loss due to agriculture. All maps are provided in the 10 km resolution projected in the geographical coordinates.

![Global map of SOCS loss due to land use](https://github.com/whrc/Soil-Carbon-Debt/blob/master/SOCS/Figure_SOC_debt_due_to_land_use.png "Global distribution of cropping and grazing in 2010 from HYDE v3.2 (above) and modeled soil organic carbon (SOC) change in the top two meters (below). Color gradients (image above) indicate proportion of grid cell occupied by given land use. Legend presented (image below) as histogram of SOC loss (Mg C / ha) with positive indicating loss and negative values depicting net gains in SOC.")

Raster files of predictions are organized by folder:

* `OCD` = Organic carbon density (kg C / cubic-m) for 0, 30, 100 and 200 cm depths predicted at various years,
* `SOCS` = Soil organic carbon stocks (Mg C / ha) for 0-30, 0-100 and 0-200 cm depths predicted at various years,
* `abs_error` = Absolute error in prediction of OCD (kg C / cubic-m) at 0, 30, 100 and 200 cm depths,

Meta-analysis data have been deposited on [Harvard University's website](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QQQM8V).

Media coverage:

* August 21, 2017 (Woods Hole Research Centre): [*New study finds soil carbon losses nearly equal to total emissions from deforestation*](http://whrc.org/new-study-finds-soil-carbon-losses-nearly-equal-to-total-emissions-from-deforestation/)
* August 22, 2017 (Thomson Reuters): [*Agriculture a culprit in global warming, says U.S. research*](http://www.reuters.com/article/us-global-climatechange-agriculture-idUSKCN1B20TR)
* August 23, 2017 (Washington Post): [*This is why when you talk about climate change, you can’t ignore agriculture*](https://www.washingtonpost.com/news/energy-environment/wp/2017/08/23/this-is-why-when-you-talk-about-climate-change-you-cant-ignore-agriculture/?utm_term=.010503c06208)
* August 25, 2017 (CarbonBrief.org): [*World’s soils have lost 133bn tonnes of carbon since the dawn of agriculture*](https://www.carbonbrief.org/worlds-soils-have-lost-133bn-tonnes-of-carbon-since-the-dawn-of-agriculture)


Please cite as:

* Sanderman, J., Hengl, T., Fiske, G., 2017 [**"The soil carbon debt of 12,000 years of human land use"**](http://dx.doi.org/10.1073/pnas.1706103114), PNAS 114(36): 9575–9580. doi:10.1073/pnas.1706103114
