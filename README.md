## Build the docker file in the folder with the dockerfile
docker build -t calc-gsva-scores .

## Run the analysis
### input variables: expression data file, geneset file, number of cores, number of bootstraps, prefix output files
sudo docker run --rm -v ~/GSVA/acc_files/:/study:ro -v ~/GSVA/results_acc_tcga/:/outdir calc-gsva-scores ./calculate_GSVA_scores.R /study/data_RNA_Seq_v2_expression_median.txt /study/meta_RNA_Seq_v2_expression_median.txt /study/msigdb.v6.0.entrez.gmt 5 500 /outdir/acc_tcga


