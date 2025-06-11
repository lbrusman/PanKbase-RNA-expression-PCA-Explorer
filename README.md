# ![](files/PanKbase_logo-black-tagline.svg) &nbsp;&nbsp; <span style = "font-family:'google', 'Open Sans', sans-serif; font-weight: 600;" > PCA Explorer
![](files/bitmap10.png)

 This is the Shiny application that runs the PCA Explorer Tool for [PanKbase](https://pankbase.org). See it live here.
 
 All data used in this app are available via the [PanKbase Integrated Cell Browser](https://dev.pankbase.org/single-cell.html?datasetId=islet_of_Langerhans_scRNA_v3-3). 
 These data are from the "Single cell expression map of pancreatic islets using data from HPAP, IIDP, and Prodo" (version v3).
 
 How to run this app yourself:
 
 - First, run [get_pseudobulk_cpm.R](app/code/get_pseudobulk_cpm.R) to generate pseudo-bulk gene expression matrices
 
 - Then, generate PCA results using [generate_PCA_results.R](app/code/generate_PCA_results.R)
 
 - Finally, [app.R](app/app.R) can be run to visualize results
