# MFAssignR
## Package Overview and References
The MFAssignR package was designed for multi-element molecular formula (MF) assignment of ultrahigh resolution mass spectrometry measurements. A number of tools for internal mass recalibration, MF assignment, signal-to-noise evaluation, and unambiguous formula selections are provided. This package contains MFAssign(), MFAssignCHO(), MFAssignAll(), SNplot(), SNcutCheck(), MFRecalList(), MFRecalCheck(), and IsoFiltR() described in the sections below. Caution with parameter settings and output evaluation is recommended.

## Molecular Formula Assignment
The molecular formula assignment algorithm in MFAssign was adapted from the low mass moiety CHOFIT assignment algorithm developed by Green and Perdue (2015). Briefly, the CHOFIT algorithm uses low mass moieties such as CH4O-1 and C4O-3 to move around in the O/C and H/C space to assign masses to CHO formulas. These low mass moieties efficiently assign CHO formulas without conventional loops. Additional combinatorial assignments with various heteroatoms are made using nested loops that subtract the mass of a heteroatom from the measured ion mass, creating a CHO “core” mass, which can then be assigned using the low mass moiety CHOFIT approach. This is further explained in Green and Perdue (2015) and Perdue and Green (2015). 

### MFAssign Functions
In total there are 3 versions of MF Assign, including MFAssign(), MFAssignCHO(), and MFAssignAll(). Where MFAssign() and MFAssignAll() include external nested loops to assign additional heteroatoms, as described in Green and Perdue (2015) and MFAssignCHO() does not. 

#### MFAssign()
Using the low mass moiety and combinatorial assignment approach, MFAssign() can be used to assign molecular formulas with 12C, 1H, and 16O and a variety of heteroatoms and isotopes, including 2H, 13C, 14N, 15N, 31P, 32S, 34S, 35Cl, and 37Cl. It can also assign Na+ adducts, which are common in positive ion mode. Due to the increasing number of chemically reasonable molecular formulas with the increasing number of possible elements and increasing molecular weight, the output will provide a list of ambiguous and unambiguous molecular formulas. 

Advanced Kendrick mass and z* sorting tools are used to reduce the number of ambiguous molecular formulas in MFAssign(). First, Kendrick mass defect (KMD) and z* values are calculated with a CH2 Kendrick base to sort the measured masses into CH2 homologous series (Stenson et al., 2003). The function then selects 1 to 3 members of each CH2 homologous series with masses below the user defined cutoff and attempts to assign molecular formulas. The ambiguous MF are then returned to the unassigned list. Then, the unambiguous MF are used as seeds for additional assignments using CH2, O, H2, H2O, and CH2O molecular formula extensions (Kujawinski and Behn, 2006). To do the formula extensions the KMD and z* values for each of these bases are calculated and then used to assign MF through the addition or subtraction of the series bases. MFAssign() (and MFAssignCHO()) tracks how many different “paths” can be used to assign each formula and if a single mass has multiple formulas, the function will choose the formula that has the largest number of paths that intercept with it. For example, if a single mass has two possible MF and one has 20 potential “paths” to it, while the other has 4, the function will choose the MF with 20 paths. Work is ongoing to track these paths and the removed MF in the data frame output of these functions. Overall, the multi-path MF extension approach greatly reduces the number of ambiguous assignments and provides an increased level of confidence in the final MF list because the MF are related to unambiguous MF assigned below the user defined cut point.

#### MFAssignCHO()
MFAssignCHO() is a simplified version of MFAssign() used to assign CHO molecular formulas only. MFAssignCHO() runs faster than MFAssign() and is best used as a preliminary formula assignment step prior to the selection of recalibrant ions in conjunction with MFRecalList() and MFRecalCheck(), which are described below. 

#### MFAssignAll()
MFAssignAll() uses the low mass moiety and combinatorial assignment approach with a simplified MF extension approach. However, only the CH2 and H2O Kendrick bases are used for MF assignment. This function results in a significantly higher number of ambiguous MF and is best used to assign previously unassigned masses after MFAssign() or on short mass lists without a complex mixture.

### Preliminary Isotope Filtering
The IsoFiltR() function can identify and pre-screen many of the 13C isotope masses, which can lower the number of peaks assigned an incorrect molecular formula. This function operates on a two column dataframe using the same structure as the MFAssign() function. IsoFiltR() works by finding the 13C Kendrick mass defects and then using those in conjunction with the CH2 z* values to differentiate the possible monoisotope/polyisotope pairs. For QA the function checks that the "polyisotopic 13C" peak is less than a certain fraction of the abundance of the "monoisotopic" peak. There are two data frame outputs from this function, "Mono" and "Iso". The "Mono" dataframe contains the masses identified as likely monoisotopes from monoisotope/polyisotope pairs and the remaining unmatched masses The "Iso" dataframe contains only the masses that were identified as polyisotopic 13C masses from the monoisotope/polyisotope pairs.
 
When the two dataframe outputs from IsoFiltR() are put into MFAssign(), the function will match the assigned monoisotopic mass peaks to their corresponding polyisotopic mass peak. Additional work is needed to use the isotopes to reduce ambiguous molecular formula assignments assigned to a single mass. IsoFiltR() should not be considered as definitive proof of the presence or absence of polyisotopic 13C molecular formulas, but it does provide some ability to identify theses masses and limit the chances that they are incorrectly assigned with a monoisotopic formula.

### Molecular Formula Quality Assurance
MFAssign() includes a number of quality assurance (QA) steps to ensure that the assigned molecular formulas are chemically reasonable. Relatively lenient default settings are provided to avoid removing chemically reasonable ambiguous molecular formula assignments. Many of these parameters are customizable, including  DBE-oxygen limits (Herzsprung et al. Anal. and Bioanal. Chem. 2014), oxygen-to-carbon ratio limits, hydrogen-to-carbon ratio limits, and minimum number of oxygen limits. The Hetcut parameter can be used to select the molecular formula with the lowest number of heteroatoms if more than one molecular formula is assigned to a single mass (Ohno and Ohno, 2013). The NMScut parameter identifies the CH4 vs O exchange series in each nominal mass as described in Koch et al. (2007), which can be used to limit ambiguous assignments. Additional non-adjustable QA parameters are used in MFAssign(), including the nitrogen rule, large atom rule, and the maximum number of hydrogen rule, maximum DBE rule (Lobodin et al., 2012), and the Senior rules (Kind et al. 2007).

## Signal-to-Noise Assessment
Signal-to-Noise level assessment can be accomplished using the SNcutCheck() and SNplot() functions, which are based on the method developed by Zhurov et al. (2014). This method uses the histogram distribution of the natural log intensities in a raw mass spectrum to determine the point where noise peaks give way to analyte signal. Additionally, the SNplot() function allows qualitative assessment of the effectiveness of the chosen S/N cut.

The SNcutCheck() function is used to estimate the signal-to-noise cut for the raw mass spectrum output mass list. The recommended signal-to-noise cut is reported in the console and a histogram showing the distribution of natural log of the abundance (or intensity) values is generated. The cut point is denoted using red for values below the cut point and blue for those above. This function should be run on the mass list prior to molecular formula assignment with MFAssign(). Setting a reasonable S/N cut point greatly increases the speed of the function and improves the output quality.
 
The SNplot function is used to show the mass spectrum with the formulas below the cut and above the cut denoted by the same colors as in the histogram plot from SNcutCheck.

## Internal Mass Recalibration
MFRecalList() and MFRecalCheck() are functions pertaining to the internal mass recalibration method adapted from Kozhinov et al. (2013) using a polynomial central moving average to estimate the weights used to recalibrate the masses. The function MFRecalList() can be used with the output of MFAssign() or MFAssignCHO() to generate a data frame containing potential recalibrant CH2 homologous series. There are a variety of metrics included in the output of this function to aid the user in picking suitable recalibrant series, these are described in greater detail in the example of MFRecalList() below. The user can select up to 10 homologous series as inputs for the mass recalibration with MFRecalCheck(). MFRecalCheck() takes the chosen series and then uses the H2 and O KMD and z* series to identify additional formulas that are related to the chosen recalibrants. To avoid recalibration problems, the function assigns the user defined number of  tallest peaks within a user defined mass range “bin” as recalibrants. For example, if the bin width is set at 20 and the number of peaks is set at 2, the function will select the two tallest peaks within each 20 m/z window across the range of the spectrum. The user can set the mass window bins and the number of peaks that will be chosen as recalibrants within each bin. This function then recalibrates the mass list and generates a plot with the input recalibrant series highlighted for a qualitative look at their overall quality.

## Recommended Order of Operations
The functions will be described in the order that they are most effectively used. The functions do not have to be run in this order, but the best results will likely be obtained in this way. A list of the functions in the recommended order is given below:
 
1. Run SNcutCheck() to determine the signal-to-noise cut point for the data
 
2. Check effectiveness of S/N cut point using SNplot()
 
3. Use IsoFiltR() to identify potential 13C isotope peaks
 
4. Using the S/N cut point determined by SNcutCheck(), and the two dataframes output from IsoFiltR(), run MFAssignCHO() to assign CHO formulas for use in identifying recalibrant ions.
 
5. Use MFRecalList() to generate a list of potential recalibrant series.
 
6. After choosing a few recalibrant series, use MFRecalCheck() to check whether they are good recalibrants and recalibrate the mass lists using those recalibrants.
 
7. Use MFAssign() with the recalibrated mass lists to assign molecular formulas to the data.
 
8. Check the output plots from MFAssign() to check the quality of the assignment.

## References
Green, N. W. and Perdue, E. M.: Fast Graphically Inspired Algorithm for Assignment of Molecular Formulae in Ultrahigh Resolution Mass Spectrometry, Anal Chem, 87(10), 5086–5094, doi:10.1021/ac504166t, 2015.

Gross, J. H.: Mass Spectrometry, doi:10.1007/978-3-319-54398-7, 2017. 

Herzsprung, P., Hertkorn, N., Tümpling, W. von, Harir, M., Friese, K. and Schmitt-Kopplin, P.: Understanding molecular formula assignment of Fourier transform ion cyclotron resonance mass spectrometry data of natural organic matter from a chemical point of view, Anal Bioanal Chem, 406(30), 7977–7987, doi:10.1007/s00216-014-8249-y, 2014.

Koch, B. P., Dittmar, T., Witt, M. and Kattner, G.: Fundamentals of Molecular Formula Assignment to Ultrahigh Resolution Mass Data of Natural Organic Matter, Anal Chem, 79(4), 1758–1763, doi:10.1021/ac061949s, 2007.

Kozhinov, A. N., Zhurov, K. O. and Tsybin, Y. O.: Iterative Method for Mass Spectra Recalibration via Empirical Estimation of the Mass Calibration Function for Fourier Transform Mass Spectrometry-Based Petroleomics, Anal Chem, 85(13), 6437–6445, doi:10.1021/ac400972y, 2013.

Kujawinski, E. B. and Behn, M. D.: Automated Analysis of Electrospray Ionization Fourier Transform Ion Cyclotron Resonance Mass Spectra of Natural Organic Matter, Anal Chem, 78(13), 4363–4373, doi:10.1021/ac0600306, 2006.

Lobodin, V. V., Marshall, A. G. and Hsu, C. S.: Compositional Space Boundaries for Organic Compounds, Anal Chem, 84(7), 3410–3416, doi:10.1021/ac300244f, 2012.

Ohno, T. and Ohno, P. E.: Influence of heteroatom pre-selection on the molecular formula assignment of soil organic matter components determined by ultrahigh resolution mass spectrometry, Anal Bioanal Chem, 405(10), 3299–3306, doi:10.1007/s00216-013-6734-3, 2013.

Perdue, E. M. and Green, N. W.: Isobaric Molecular Formulae of C, H, and O: A View from the Negative Quadrants of van Krevelen Space, Anal Chem, 87(10), 5079–5085, doi:10.1021/ac504165k, 2015.

Zhurov, K. O., Kozhinov, A. N., Fornelli, L. and Tsybin, Y. O.: Distinguishing Analyte from Noise Components in Mass Spectra of Complex Samples: Where to Cut the Noise, Anal Chem, 86(7), 3308–3316, doi:10.1021/ac403278t, 2014. 
