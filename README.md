# scRNA-Seq
# Author: Brian Wray

### General instructions

- Open a window and navigate to a quest analytics node
- File > New Project... 
- (In new project wizard) Version Control > Git 
- Repository URL: https://github.com/brianFSM/scRNA-seq_seurat_V5.git

the renv.lock that is downloaded

Go back to your project in the analytics node, and at the console:

renv::init(bare=TRUE)
renv::restore()

If you've ever installed the libraries before, they will link from the cache in a flash. Otherwise they will be downloaded and compiled. 

If this is the first time you've ever run these reports, you will need to run the following from the console in the analytics node in order to get pdfs to render (the reports are now pdfs, not html files):

> install.packages("tinytex")
> tinytex::install_tinytex()   # downloads + sets up ~/.TinyTeX
> tinytex::is_tinytex()        # should return TRUE

Finally, fill out your config file. 

You can now run multiple templates in succession. 
To run the templates from the command line on quest, do something like this:

$ sbatch run_templates.sh scRNA_part1_QC.Rmd scRNA_part2a_integration.Rmd

In this example, the output will be scRNA_part1_QC.pdf and scRNA_part2a_integration.pdf

