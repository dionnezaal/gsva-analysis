## Building the docker file in folder where the docker file is located
``` docker build -t calc-repr-gsva-scores . ```

## Running the docker to calculate representative score for a study
``` sudo docker run --rm -v $GSVA_FILE_DIR:/study:ro -v $OUTDIR:/out calc_repr_scores ./calc_representative_geneset_score.R /study/$GSVA_SCORE_FILE /out/$PREFIX_OUTPUT ```
