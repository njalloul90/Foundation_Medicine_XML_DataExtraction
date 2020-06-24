# Foundation_Medicine_Data_Extraction
R Shiny application to parse and extract relevant information from Foundation Medicine tumor sequencing xml files. 

# Foundation Medicine's FoundationOne CDx test

The FoundationOne CDx test is a companion diagnostic test for solid tumors. The test is designed to provide clinically actionable information from solid tumor sequencing. The full description of the test, including the gene list, biomarker description, and technical information can be found at:

https://www.foundationmedicine.com/genomic-testing/foundation-one-cdx

# Descripton of Foundation Medicine's tumor sequencing xml files

A detailed XML (Extensible Markup Language) data file is usually provided by Foundation Medicine (or can be requested) along with a pdf report that presents the major findings (general infromation, diagnosis, identified genomic alterations, biomarkers, therapeutic interventions and available clinical trials associated with the identified alterations, and variants of unknown significance). A sample of the pdf report canbe found at:

https://assets.ctfassets.net/vhribv12lmne/P1UbtVjOoeAcaOCWoWQkW/7e3d66f7396466156e8eb4b27d0d471b/F1CDx_SampleReport.pdf

Most bioinformatic analysis as well as clinical investigations into the patient's tumor profile require further information that is left out of the pdf report delivered to physicians. These include the variant allele frequencies, sequencing depth for each variant, chromosomal loci of the mutation in order to determine its effect on the gene... 
For that purpose, this tool is presented to to provide clinicians with a simple user interface able to extract all relevant information from the xml files for analysis.

# Running FoundationMedicineDataExtraction application

The application can be directly accessed through the url:

https://njalloul.shinyapps.io/FoundationMedicineDataExtraction/?_ga=2.21770171.1421588860.1593013867-628062037.1593013867

The tumor sequencing results presented in the xml files are organised in separate tags. The table below indicates the tags used to extract information in this application as well as the description of each tag:

| Tag | Description |
| --- | --- |
| ReportId | Unique report ID |
| SubmittedDiagnosis | Resulting diagnosis from the solid tumor sequencing |
| Gender | Patient's listed gender |
| DOB | Patient's date of birth |
| SpecSite or SpecFormat | Specimen site delivered for sequencing |
| CollDate | Date of sample collection at medical facility |
| ReceivedDate | Date of sample recieved by FoundationMedicine |
| CountryOfOrigin | Patient's country of origin |
| percent-tumor-nuclei | Histological tumor purity estimate |
| purity-assessment | Computational tumor purity estimate |
| tumor-mutation-burden score | Tumor mutational burden score and status |
| microsatellite-instability status | Microsatellite instability status |
| short-variants | Identified genomic alterations; attributes include allele frequencies, sequencing depth,protein change, chromosomal position and strand, functional effect, and mutation |
| copy-number-alterations | Identified copy number alterations; alterations found here can overlap with the identified genomic alterations |
| VariantProperties | Variants of unknown significance (VUS) |
| rearrangements | Identified rearrangemets; attributes include targeted gene, description, second gene, and corresponding positions |




