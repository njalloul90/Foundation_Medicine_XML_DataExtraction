# Foundation_Medicine_Data_Extraction
R Shiny application to parse and extract relevant information from Foundation Medicine tumor sequencing xml files. 

# Foundation Medicine's FoundationOne CDx test

The FoundationOne CDx test is a companion diagnostic test for solid tumors. The test is designed to provide clinically actionable information from solid tumor sequencing. The full description of the test, including the gene list, biomarker description, and technical information can be found at:

https://www.foundationmedicine.com/genomic-testing/foundation-one-cdx

# Descripton of Foundation Medicine's tumor sequencing xml files

A detailed xml file is usually provided by Foundation along with a pdf report that present the major findings (general infromation, diagnosis, identified genomic alterations, biomarkers, therapeutic interventions associated with the identified alterations, and variants of unknown significance). A sample of the pdf report canbe found at:
https://assets.ctfassets.net/vhribv12lmne/P1UbtVjOoeAcaOCWoWQkW/7e3d66f7396466156e8eb4b27d0d471b/F1CDx_SampleReport.pdf

Most bioinformatic analysis as well as clinical investigations into the patient's tumor profile require further information that is left out of the pdf report delivered to physicians. These include the variant allele frequencies, sequencing depth for each variant, chromosomal loci of the mutation in order to determine its effect on the gene... 
For that purpose, this tool is presented to to provide clinicians with a simple user interface able to extract all relevant information from the xml files for analysis.

# Running FoundationMedicineDataExtraction application

The application can be directly accessed through the url:

https://njalloul.shinyapps.io/FoundationMedicineDataExtraction/?_ga=2.21770171.1421588860.1593013867-628062037.1593013867

A sample xml file is included to test out the application.
