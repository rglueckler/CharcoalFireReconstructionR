# CharcoalAnalysis in R (v.1.0)
R script for the analysis of macroscopic charcoal records from sedimentary archives. It combines previous work (please see references listed below) to provide a quick way of evaluating fire signals within the record using both a "classic" and a "robust" approach, featuring:
- Interpolation to equal time steps
- Calculation of charcoal accumulation rates (CHAR)
- Separation into background and peak components
- Signal-to-noise index (SNI) for the peak component
- Identification of fire episodes and their SNI-derived quality
- Robust charcoal analysis approach, taking into account age and proxy uncertainties
- Automated basic output plots (.pdf)
- Automated output information (.txt) about settings and results, including fire return intervals (FRI)

This script was created in R v.4.0.2 and utilizes packages "paleofire", "locfit" and "mixtools". No warranty, please use at your own risk.

# How to use the script
After downloading the and extracting the files to a working directory, the "CharcoalAnalysis_ExampleRun.R" file provides a quick way to test the script with an exemplary charcoal record created with random numbers. The according file, which might also be used as a template for custom data, is located in: "/Records/RandomCore/RandomCore.csv". For easy handling within the script new records should be added in folders sharing the record's name, e.g. "/Records/NewRecord1/NewRecord1.csv". After starting the "CharcoalAnalysis_ExampleRun.R" script, follow the steps outlined with #comments. After declaring the name of the record (e.g. record_name = "NewRecord1") and looking through the other input variables that could be adjusted, the script could be run in its entirety. Output files will be saved within the respective folder of a record. For more details on the methods and functions used, please refer to references mentioned below. 
Custom data for each sample of a record should be in a .csv table with the following columns: 1: "depth_top" (depth at top/upper end of sample interval), 2: "depth_bot" (depth at bottom/lower end of sample interval), 3: "age_top" (age at top/upper end of sample interval), 4: "age_bot" (age at bottom/lower end of sample interval), 5: "age_uncert" (age uncertainty, e.g. 1 sig range of estimated ages), 6: "vol" (volume, used in cm³), 7: "proxy" (i.e. absolute charcoal counts), 8: "proxy_uncert" (proxy uncertainty or error, i.e. deviation estimated by counts from duplicate samples). 

# References
For the original version of "CharAnalysis", please check the respective repository by Philip Higuera (https://github.com/phiguera/CharAnalysis) and these publications:

<blockquote>Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, T.A. Brown (2009): Vegetation mediated the impacts of postglacial climatic change on fire regimes in the south-central Brooks Range, Alaska. Ecological Monographs 79: 201-219</blockquote>

<blockquote>Higuera, P.E., D.G. Gavin, P.J. Bartlein and D.J. Hallett (2010): Peak detection in sediment-charcoal records: impacts of alternative data analysis methods on fire-history interpretations. International Journal of Wildland Fire 19: 996-1014</blockquote>

For the translation of "CharAnalysis" into R used in the present script, featuring a detailed description of the robust analysis approach, please see:

<blockquote>Dietze, E., D. Brykała, L.T. Schreuder, K. Jażdżewski, O. Blarquez, A. Brauer, M. Dietze, M. Obremska, F. Ott, A. Pieńczewska, S. Schouten, E.C. Hopmans, M. Słowiński (2019): Human-induced fire regime shifts during 19th century industrialization: A robust fire regime reconstruction using northern Polish lake sediments. PloS one 14(9)</blockquote>

For further detail and the script on the SNI, please see:

<blockquote>Kelly, R., P.E. Higuera, C. Barrett, F. Hu (2011): A signal-to-noise index to quantify the potential for peak detection in sediment-charcoal records. Quaternary Research 75(1): 11-17</blockquote>

# Citation
If you use this script for a publication, please cite the literature mentioned above.

# Contact
For questions or comments/feedback, feel free to contact ramesh.glueckler@awi.de
