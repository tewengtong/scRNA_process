# scRNA_process
Here includes raw data process, annotation, and how to further analysis scRNA data for sQTL

## CellRanger

At first, you need to download reference genome and its gtf file, since our data comes from human, you can download them from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

Next, you can use the count function to generate the files for Seurat:

```
nohup bash cellranger_count.sh &
```

## Annotation

Here, we recommend you to annotate the cell cluster manually, and this R script includes how to process data coming from cellranger and get the files for annotation.
