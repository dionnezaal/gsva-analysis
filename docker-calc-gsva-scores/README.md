## Build the docker file in the folder with the dockerfile
``` docker build -t calc-gsva-scores . ```

## Run the analysis
input variables: expression data file, geneset file, number of cores, number of bootstraps, prefix output files

``` sudo docker run --rm -v ~/$STUDY_DIRECTORY/:/study:ro -v ~/$OUTPUT_DIRECTORY/:/outdir calc-gsva-scores ./calculate_GSVA_scores.R /study/$EXPRESSION_FILE /study/$META_EXPRESSION_FILE /study/$GENESET_FILE $NUMBER_OF_CORES $NUMBER_OF_BOOTSTRAPS /outdir/$PREFIX_OUTPUT_FILES ```


