## activate snakemake
conda activate snakemake
mamba activate snakemake

## buld docker container
sudo docker build -t vibaotram/simulation:1.0 .
sudo docker run -it vibaotram/simulation:1.0
sudo docker push vibaotram/simulation:1.0

## run the pipeline
snakemake --sdm apptainer conda --configfile config.yaml