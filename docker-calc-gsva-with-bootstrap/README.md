## Build the docker file in the folder with the dockerfile
``` docker build -t calc-gsva-bootstrap . ```

## Run the analysis
input variables: expression data file, geneset file, meta expression file, prefix output files, number of cores, number of bootstraps

``` sudo docker run --rm -v ~/$STUDY_DIRECTORY/:/study:ro -v ~/$OUTPUT_DIRECTORY/:/outdir calc-gsva-bootstrap ./calc_GSVA_with_bootstrap_function.R /study/$EXPRESSION_FILE /study/$GENESET_FILE /study/$META_EXPRESSION_FILE /outdir/$PREFIX_OUTPUT_FILES $NUMBER_OF_CORES $NUMBER_OF_BOOTSTRAPS ```

## The bootstrapping method is explained in: GSVA-bootstrap-method.md
