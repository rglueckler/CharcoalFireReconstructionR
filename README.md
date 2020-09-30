# CharcoalAnalysis in R
R script for the analysis of macroscopic charcoal records from sedimentary archives. It combines previous work (please check publications listed below) to provide a quick way of evaluating fire signals within the record, featuring:
- interpolation to equal time steps
- calculation of charcoal accumulation rates (CHAR)
- determination of background and peak components
- a signal to noise-ratio index (SNI) for the peak component
- identification of fire episodes and their SNI-derived quality
- an additional robust analysis approach taking into account age and proxy uncertainties
- automated output plots (.pdf)
- automated output information (.txt) about settings and results, including fire return intervals (FRI)

This script was created in R v.4.0.2 and utilizes packages "paleofire", "locfit", "mixtools" and "stats". Use without warranty.

# How to use the script
Download the .zip file and unpack into a new folder. Open the R file with the R tool of your choice. Your data should be in a .csv table with the following columns: 
1: depth (top), 2: depth (bottom), 3: age (top), 4: age (bottom), 5: volume of sediment, 6: number of counted charcoal particles.

# Literature
For the original version of "CharAnalysis", please check the respective repository by Philip Higuera (https://github.com/phiguera/CharAnalysis) and these publications:

<blockquote>Higuera, P.E., L.B. Brubaker, P.M. Anderson, F.S. Hu, T.A. Brown (2009): Vegetation mediated the impacts of postglacial climatic change on fire regimes in the south-central Brooks Range, Alaska. Ecological Monographs 79: 201-219</blockquote>

<blockquote>Higuera, P.E., D.G. Gavin, P.J. Bartlein and D.J. Hallett (2010): Peak detection in sediment-charcoal records: impacts of alternative data analysis methods on fire-history interpretations. International Journal of Wildland Fire 19: 996-1014</blockquote>

For the translation of "CharAnalysis" into R used in the present script, featuring a detailed description of the robust analysis approach, please see:

<blockquote>Dietze, E., (2019):</blockquote>

For further detail and the script on the SNI, please see:

<blockquote>Kelly X., (2011):</blockquote>

# Citation
If you use this script for a publication, please cite the literature mentioned above.

# Contact
For questions or comments/feedback, feel free to contact ramesh.glueckler@awi.de
