# MFAssignR


##Package Overview and References

The MFAssignR package was designed for multi-element molecular formula (MF) assignment of ultrahigh resolution mass spectrometry measurements. A number of tools for internal mass recalibration, MF assignment, signal-to-noise evaluation, and unambiguous MF assignments are provided. This package contains MFAssign(), MFAssign_RMD(), MFAssignCHO(), MFAssignCHO_RMD(), SNplot(), HistNoise(), KMDNoise(), RecalList(), Recal(), and IsoFiltR() described in the sections below. Note, the functions with “RMD” were designed to be run within an R Markdown file and are otherwise identical to the corresponding non-”RMD” versions. To learn more, please see the section titled “Semi-Automated MFAssignR Functions” in the User Manual.  User caution with the function parameter settings and output evaluation is required; thus, several function outputs are provided to assist the user with these evaluations.

##Molecular Formula (MF) Assignment

The MF assignment algorithm in MFAssign was adapted from the low mass moiety CHOFIT assignment algorithm developed by Green and Perdue (2015). Briefly, the CHOFIT algorithm uses low mass moieties such as CH4O-1 and C4O-3 to move around in the O/C and H/C space to assign MF with C, H, and O without conventional loops. The MFAssignCHO() function uses the CHOFIT strategy to assign MF with C, H, and O. Additional combinatorial assignments with various heteroatoms are made using nested loops that subtract the mass of a heteroatom from the measured ion mass, creating a CHO “core” mass, which can then be assigned using the low mass moiety CHOFIT approach. The MFAssign() function uses this latter approach with several additional heteroatoms. Further information is available in Green and Perdue (2015) and Perdue and Green (2015). 

###MFAssign()

Using the low mass moiety and combinatorial assignment approach, MFAssign() can be used to assign MF with 12C, 1H, and 16O and a variety of heteroatoms and isotopes, including 2H, 13C, 14N, 15N, 31P, 32S, 34S, 35Cl, 37Cl,and 19F. It can also assign Na+ adducts, which are common in positive ion mode. Due to the increasing number of chemically reasonable MF with the increasing number of possible elements and increasing molecular weight, the output will provide a list of ambiguous and unambiguous MF. 

In MFAssignR, we use a ‘de novo’ concept for MF assignment, where ‘de novo’ means the first in series. This approach takes advantage of the naturally occurring mass spectral patterns typically observed in natural organic matter. The most frequent nominal mass difference patterns include: 2, 14, and 16 that correspond to H2, CH2, and O. Thus, these patterns are used to restrain the number of chemically reasonable MF assigned to masses above the user defined ‘de novo’ cutoff (e.g., 300). In MFAssign(), this is done using Kendrick mass defects and z* sorting. First, Kendrick mass defects (KMD) and z* values are calculated with a CH2 Kendrick base to sort the measured masses into CH2 homologous series (Stenson et al., 2003). The function then selects 1 to 3 members of each CH2 homologous series with masses below the user defined cutoff and attempts to assign MF. The ambiguous MF are then returned to the unassigned list. Then, the unambiguous MF are used as seeds for additional assignments using CH2, O, H2, H2O, and CH2O MF extensions (Kujawinski and Behn, 2006). To do the formula extensions, the KMD and z* values for each of these bases are calculated and then used to assign MF through the addition or subtraction of the series bases. 

MFAssign tracks how many “paths” can be used to assign each MF and if a single mass has multiple MF. By default, the function will choose the MF that has the largest number of paths that intercept with it. For example, if a single mass has two possible MF and one has 20 potential “paths” to it, while the other has 4, the function will choose the MF with 20 paths. Work is ongoing to track these paths and the associated MF in the data frame output of these functions. Overall, the multi-path MF extension approach greatly reduces the number of ambiguous assignments and provides an increased level of confidence in the final MF list because the MF are related to unambiguous MF assigned below the user defined cutoff. To reduce the number of ambiguous sulfur assignments, sulfur containing MF used as seeds must be unambiguous and have a matching 34S peak. 

To allow ambiguity in the formula assignments there is the "Ambig" parameter which can be turned "on" or "off". This option turns off the path frequency prioritization step for the formula assignments as described above, which allows all chemically reasonably MF assignments to be retained for each mass. Additionally, an "MSMS" parameter is available, which can be used to assign MF in a data set that is not very continuous (e.g., MS/MS data). In this case, no pre-filtering of the masses below the “DeNovo” threshold is done, meaning that all masses below the threshold will be assigned directly. This causes the function to run somewhat slower, but can improve assignment coverage in some situations. These parameters replace the MFAssignAll() and MFAssignMSMS() functions from previous versions (<= v.0.0.3).

###MFAssignCHO() 

MFAssignCHO() is a simplified version of MFAssign() used only to assign MF with CHO elements. MFAssignCHO() runs faster than MFAssign() and can be used for preliminary MF assignments prior to the selection of internal recalibration ions in conjunction with RecalList() and Recal(), which are described below. 

##Isotope Filtering

The IsoFiltR() function can identify prospective 13C and 34S isotope masses. This is done to avoid incorrect monoisotopic MF assignments. This function operates on a two column data frame using the same structure as the MFAssign() function. 

IsoFiltR() identifies potential isotope masses using a four-step identification method. 

1. First the mass list is transformed to identify mass difference pairs appropriate for the element under investigation (delta mass for C (1.003355) or S (1.995797), with +/- 5 ppm mass error). Only those that meet this criterion move on to step 2.

2. Using the mass difference between 12C/13C (1.003355) or 32S/34S (1.995797), the KMD value can be calculated for a specific isotope. This means that the 12C (32S) monoisotopic peak will be in a KMD homologous series with its matching 13C (34S) isotopic peak, analogous to homologous series of CH2. If the KMD values are equivalent for the candidate pair, the masses can be considered to be in a series and the pair will move on to the third step. The equations for 13C are: KM = 1/1.003355 * m/z and KMD = nominal mass - KM. Then, 2/1.995797 replaces 1/1.003355 for 34S.

3. Isotope pairs are separated using a “Resolution Enhanced KMD” approach adapted from Zheng et al. 2019. Resolution enhanced KMD values are calculated by dividing the mass of some homologous series base (in this case CH2) by an integer that is experimentally determined to accomplish the desired separation. This value is then used in the typical KM and KMD calculation in order to calculate the “resolution enhanced” KMD (re-KMD). For example, the integer 21 is used to adjust the CH2 base mass in the following KMD calculation: BaseMass_adj = 14.01565 / 21  and then re-KM = (round(BaseMass_adj) / BaseMass_adj) * m/z, followed by re-KMD = round(re-KM) - re-KM yields a resolution enhanced KMD. 

For 13C, the integer 21 is used in the resolution-enhanced KMD, while for 34S it is 12. Then, the masses that are 12/13 C or 32/34 S pairs will have specific re-KMD difference values, which are used to select the pairs of masses that are most likely to be isotope pairs. The re-KMD differences (polyisotope – monoisotope) are both positive and negative because the re-KM and re-KMr values were rounded off. The values are -0.291 and 0.709 for 32/34 S and -0.496 and 0.503 for 12/13 C. If the masses meet these criteria, they can move on to step four. Using CH2 KMD values that are divided by an experimentally derived integer, the isotope pairs are separated into two specific values. If the difference in the enhanced KMD for the candidate pair matches one of those values, it will move to the fourth step.

4. The abundance ratios are used to constrain the remaining isotope pairs to ensure that the isotope masses are not too large or too small relative to the intensity of the monoisotopic peak. The limits on this are loose due to the variation in the polyisotope abundance with analyte signal (similar to isotope dilution) as observed in ultrahigh resolution Orbitrap and FT-ICR measurements.

The candidate pairs that make it through these four steps are put into two data frames, Mono and Iso, which contain the monoisotopic and isotopic masses respectively. Then all of the masses that were not flagged as possible mono/iso pairs are returned to the Mono output data frame. In complex mixtures, some masses can be flagged as both monoisotopic and isotopic. In these cases, the masses are included in both outputs and are classified as either monoistopic or isotopic after the MF assignment.
  
When the two data frame outputs from IsoFiltR() are put into MFAssign(), the function will match the assigned monoisotopic masses to their corresponding isotopic masses. Additional work would be needed to use the isotopes to reduce ambiguous MF assignments assigned to a single mass. Thus IsoFiltR() should not be considered as definitive proof of the presence or absence of 13C or 34S in a MF, but it does identify most MF with these naturally occurring isotopes and limit the chances that they are incorrectly assigned with a monoisotopic MF.

##Molecular Formula (MF) Quality Assurance

MFAssign() includes a number of quality assurance (QA) steps to ensure output of chemically reasonable MF. In general, the default settings are relatively lenient to yield a wide range of chemically reasonable MF for a broad range of experiments. However, many of the parameters are customizable, including DBE-O limits (Herzsprung et al. Anal. and Bioanal. Chem. 2014), O/C ratio limits, H/C ratio limits, and minimum number of O. The HetCut parameter can be used to select the MF with the lowest number of heteroatoms, if more than one MF is assigned to a single mass (Ohno and Ohno, 2013). The NMScut parameter identifies the CH4 vs O exchange series in each nominal mass as described in Koch et al. (2007), which can be used to limit ambiguous assignments. Additional non-adjustable QA parameters are used in all of the MFAssign functions, including the nitrogen rule, large atom rule, the maximum number of H rule, maximum DBE rule (Lobodin et al., 2012), and the Senior rules (Kind et al. 2007).

##Noise Assessment

Noise level assessment can be accomplished using the either the HistNoise() or KMDNoise() functions in conjunction with the SNplot() functions. The HistNoise() method is based on the method developed by Zhurov et al. (2014), and KMDNoise() is a new custom method based on our observations of raw data Kendrick mass defect analysis. 

The Zhurov et al. (2014) method uses a histogram distribution of the natural log intensities in the measured raw mass spectrum to determine the point where noise peaks give way to analyte signal. The HistNoise() function attempts to identify this point and reports the noise level so that the signal-to-noise threshold can be determined. The threshold is shown in the output plot with red and blue colors, where red indicates noise. If the function does not predict a reasonable noise level, the threshold can be set manually by the user. We frequently observed this function to fail to separate the distributions when the analyte signal tapers into the noise. For this reason, we developed the KMDNoise() function described below.

The KMDNoise() method is based on the observation that the CH2 based KMD values of noise  and analyte masses are naturally separated in a KMD plot, allowing the function to select a region with only noise to calculate the average intensity. We refer to this as the KMD slice method. In principle, this is similar to what was described in Reidel and Dittmar (2014), but instead of using a static range of normal mass defects (0.3-0.9), our method uses a mass dependent KMD region, which avoids potentially doubly charged masses with a mass defect of ~0.5, which would be considered as noise in the Reidel and Dittmar method.  Additionally, the user can set limits on the mass range to use to estimate the noise, if that is necessary to avoid specific high intensity peaks.

At least one of these noise estimation functions should be run on the mass list prior to MF assignment with MFAssign() or isotope filtering with IsoFiltR(). Setting a reasonable S/N threshold greatly increases the speed of the functions and improves the output quality.  

The SNplot() function is used to show the mass spectrum with the masses below and above the threshold shown in the output plot with red to blue colors, where red indicates noise.

##Internal Mass Recalibration

The internal mass recalibration method in MFAssignR was adapted from Kozhinov et al. (2013) and Savory et al. (2011). It uses a polynomial central moving average to estimate the weights used to recalibrate the masses (Kozhinov et al., 2013) applied to spectral segments (Savory et al., 2011) to perform a walking recalibration. 

First, the function RecalList() can be used with the output of MFAssign() or MFAssignCHO() to generate a data frame containing potential recalibrant CH2 homologous series. There are a variety of metrics included in the output of this function to aid the user in picking suitable recalibrant mass series. Some of the more useful parameters for ensuring complete coverage of the spectrum are “Number Observed”, “Mass Range”, and “Tall Peak”. The quality of the series with regard to whether the series has “holes” in it, how tall the tallest peak is relative to other peaks in the region, and how close the two tallest peaks are to each other are estimated with the “Series Score”, “Peak Distance”, and “Abundance Score” respectively. Please see the User Manual/Vignette for more information about each parameter in RecalList(). Combined, these series should cover the full mass spectral range to provide the best overall recalibration. The best series to choose are generally long and combined have a “Tall Peak” at least every 100 m/z. 

Up to ten of these series can be chosen to be used in the Recal() function, which recalibrates the spectrum. Choosing appropriate recalibrants is a critical aspect of recalibrating a mass spectrum effectively. After selecting the recalibrant series and entering them to the Recal() function, the parameter in Recal() most likely to be changed is “mzRange” which sets the recalibration segment length and has a default value of 30. If this value does not work a warning will be printed to the R console telling the user to increase the value. Formula extension via H2 and O homologous series uses the user defined recalibrant series as a base to find addtional recalibrants. It is limited to a user defined number of steps (+/- H2 or O) and generates a pool of potential recalibrants. Formula extension occurs between the assigned unambiguous molecular formulas and then each of the potential recalibrants are also check for a matching 13C peak. If there is a matching 13C then it is also added to the pool of recalibrants to be used. This pool of recalibrants are separated into each user defined segment and used to calculate a mass error correction term based on Kozhinov et al. (2013). These mass correction terms are then used to recalibrate each segment independently, removing the systematic error present in a mass spectrum.

#Function Examples

##Recommended Order of Operations

The functions will be described in the order that they are most effectively used. The functions do not have to be run in this order, but often the best results will likely be obtained in this way. A list of the functions in the recommended order is given below:

1. Run KMDNoise() to determine the noise level for the data.
 
2. Check effectiveness of S/N threshold using SNplot().
 
3. Use IsoFiltR() to identify potential 13C and 34S isotope masses.
 
4. Using the S/N threshold, and the two data frames output from IsoFiltR(), run MFAssignCHO() to assign MF with C, H, and O to assess the mass accuracy.
 
5. Use RecalList() to generate a list of the potential recalibrant series.
 
6. After choosing recalibrant series, use Recal() to recalibrate the mass lists.
 
7. Assign MF to the recalibrated mass list using MFAssign().
 
8. Check the output plots from MFAssign() to evaluate the quality of the assignments.

The functions in the MFAssignR package were developed by adapting methods and algorithms from the peer reviewed literature. The following references were referred to in this document:

##References

Green, N. W. and Perdue, E. M.: Fast Graphically Inspired Algorithm for Assignment of Molecular Formulae in Ultrahigh Resolution Mass Spectrometry, Anal Chem, 87(10), 5086–5094, doi:10.1021/ac504166t, 2015.

Gross, J. H.: Mass Spectrometry, doi:10.1007/978-3-319-54398-7, 2017. 

Herzsprung, P., Hertkorn, N., Tumpling, W. von, Harir, M., Friese, K. and Schmitt-Kopplin, P.: Understanding molecular formula assignment of Fourier transform ion cyclotron resonance mass spectrometry data of natural organic matter from a chemical point of view, Anal Bioanal Chem, 406(30), 7977–7987, doi:10.1007/s00216-014-8249-y, 2014.

Koch, B. P., Dittmar, T., Witt, M. and Kattner, G.: Fundamentals of Molecular Formula Assignment to Ultrahigh Resolution Mass Data of Natural Organic Matter, Anal Chem, 79(4), 1758–1763, doi:10.1021/ac061949s , 2007.

Kozhinov, A. N., Zhurov, K. O. and Tsybin, Y. O.: Iterative Method for Mass Spectra Recalibration via Empirical Estimation of the Mass Calibration Function for Fourier Transform Mass Spectrometry-Based Petroleomics, Anal Chem, 85(13), 6437–6445, doi:10.1021/ac400972y, 2013.

Kujawinski, E. B. and Behn, M. D.: Automated Analysis of Electrospray Ionization Fourier Transform Ion Cyclotron Resonance Mass Spectra of Natural Organic Matter, Anal Chem, 78(13), 4363–4373, doi:10.1021/ac0600306 , 2006.

Lobodin, V. V., Marshall, A. G. and Hsu, C. S.: Compositional Space Boundaries for Organic Compounds, Anal Chem, 84(7), 3410–3416, doi:10.1021/ac300244f, 2012.

Ohno, T. and Ohno, P. E.: Influence of heteroatom pre-selection on the molecular formula assignment of soil organic matter components determined by ultrahigh resolution mass spectrometry, Anal Bioanal Chem, 405(10), 3299–3306, doi:10.1007/s00216-013-6734-3, 2013.

Perdue, E. M. and Green, N. W.: Isobaric Molecular Formulae of C, H, and O: A View from the Negative Quadrants of van Krevelen Space, Anal Chem, 87(10), 5079–5085, doi:10.1021/ac504165k, 2015.

Savory, J. J., Kaiser, N. K., McKenna, A. M., Xian, F., Blakney, G. T., Rodgers, R. P., Hendrickson, C. L., and Marshall, A. G.: Parts-Per-Billion Fourier Transform Ion Cyclotron Resonance Mass Measurement Accuracy with a "Walking" Calibration Equation, Anal Chem, 83, 1732-1736, doi:10.1021/ac102943z, 2011.

Zheng, Q., Morimoto, M., Sato, H. and Fouquet, T.: Resolution-enhanced Kendrick mass defect plots for the data processing of mass spectra from wood and coal hydrothermal extracts, Fuel, 235, 944–953, doi:10.1016/j.fuel.2018.08.085, 2019.

Zhurov, K. O., Kozhinov, A. N., Fornelli, L. and Tsybin, Y. O.: Distinguishing Analyte from Noise Components in Mass Spectra of Complex Samples: Where to Cut the Noise, Anal Chem, 86(7), 3308–3316, doi:10.1021/ac403278t, 2014. 
