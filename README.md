
# COVID19-RAGE

**scRNA-seq analysis of PBMC from COVID-19 patients of the Padova cohort**

The code and the data in this repository enables the reproduction of the analyses and plots for 
the single-cell RNA-seq data in the paper "RAGE engagement by SARS-CoV-2 enables monocyte infection and
underlies COVID-19 severity" (https://doi.org/10.1101/2022.05.22.492693).


An html version of the notebooks is accessible [here](https://GiuseppeTestaLab.github.io/covid19-RAGE/).




## Filtering, dimensionality reduction, cell annotation

Links: [jupyter notebook](01_gex_filtering_annotation.ipynb) and [html file](https://GiuseppeTestaLab.github.io/covid19-RAGE/01_gex_filtering_annotation.html).

Notebook that containes the preliminary steps of the single-cell data analysis

1) initial quality filters 
2) log-normalization of the counts and HVG selection
3) dimensionality reduction (PCA and UMAP) with Harmony integration 
4) Leiden clustering and celltype annotation
5) preliminary analysis and visualization of the annotated data




## Differential abundance analysis

Links: [html file](https://GiuseppeTestaLab.github.io/covid19-RAGE/02_abundance_analysis_cell_families.html).

RMarkdown code for the computation of the differential abundance of the
cell families that have been defined in the dataset.




## Data exploration and ingestion of external annotations

Links: [jupyter notebook](03_dataset_exploration_and_ingestion.ipynb) and [html file](https://GiuseppeTestaLab.github.io/covid19-RAGE/03_dataset_exploration_and_ingestion.html).

Plots of the relevant cell metadata and
comparison with ingested cell annotations based on
[Wilk et al.](https://www.nature.com/articles/s41591-020-0944-y)




## Pseudobulk data computation

Links: [jupyter notebook](04_pseudobulk.ipynb) and [html file](https://GiuseppeTestaLab.github.io/covid19-RAGE/04_pseudobulk.html).

Computation of the pseudobulk data for each patient sample




## Differential expression analysis

Links: [html file](https://GiuseppeTestaLab.github.io/covid19-RAGE/05_edgeR_DE.html).

Differential expression analysis of the pseudobulk data
with edgeR GL model




## RAGE pathway exploration

Links: [jupyter notebook](06_RAGE_pathway_exploration.ipynb) and [html file](https://GiuseppeTestaLab.github.io/covid19-RAGE/06_RAGE_pathway_exploration.html).

Supervised exploration of the expression of the 
genes included in the RAGE binding gene list




## External dataset summary

Links: [jupyter notebook](07_external_dataset_summary.ipynb) and [html file](https://GiuseppeTestaLab.github.io/covid19-RAGE/07_external_dataset_summary.html).

Summary of the results of the RAGE exploration
performed on the publicly available datasets of 
single-cell RNA-seq of COVID-19 patients

Detailed analisys workflows are described in the following notebooks:

- [Arunachalam et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_arunachalam_2020_processed.html)
- [Silvin et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_silvin_2020_processed.html)
- [Yao et al. 2021](https://giuseppetestalab.github.io/covid19-RAGE/08_yao_2021_processed.html)
- [Bost et al. 2021](https://giuseppetestalab.github.io/covid19-RAGE/08_bost_2021_processed.html)
- [Stephenson et al. 2021](https://giuseppetestalab.github.io/covid19-RAGE/08_stephenson_2021_processed.html)
- [Yu et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_yu_2020_processed.html)
- [Combes et al. 2021](https://giuseppetestalab.github.io/covid19-RAGE/08_combes_2021_processed.html)
- [Su et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_su_2020_processed.html)
- [Zhu et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_zhu_2020_processed.html)
- [Lee et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_lee_2020_processed.html)
- [Wen et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_wen_2020_processed.html)
- [Schulte-Schrepping et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_schulte-schrepping_2020_processed.html)
- [Wilk et al. 2020](https://giuseppetestalab.github.io/covid19-RAGE/08_wilk_2020_processed.html)



## Drug repurposing analysis

Links: [folder](https://github.com/GiuseppeTestaLab/covid19-RAGE/tree/main/10_drug_repurposing_analysis).

Analysis of the transcriptomic signatures based on the Connectivity Map (CMap) database.
This resource has been used to:

- detect whether current COVID-19 therapeutics are were able to mimic 
  the transcriptional changes occurring from in the subjects of 
  the study from admission to recovery ([01_known_drugs.R](https://github.com/GiuseppeTestaLab/covid19-RAGE/blob/main/10_drug_repurposing_analysis/01_known_drugs.R));
- determine which drugs included in the cMap database are able to 
  interfere specifically with the RAGE pathway determined from 
  Ingenuity Pathway Analysis ([02_RAGE_IPA.R](https://github.com/GiuseppeTestaLab/covid19-RAGE/blob/main/10_drug_repurposing_analysis/02_RAGE_IPA.R));
- carry our an unbiased analysis of the CMap database, seeking compounds 
  potentially capable of reverting severity signatures 
  across cell families ([03_cMAP_enrich_trt_unbiased.R](https://github.com/GiuseppeTestaLab/covid19-RAGE/blob/main/10_drug_repurposing_analysis/03_cMAP_enrich_trt_unbiased.R));
- retrieve gene expression changes upon a specific treatment 
  from level 4 cMAP data ([04_zscore_fromGSE.R](https://github.com/GiuseppeTestaLab/covid19-RAGE/blob/main/10_drug_repurposing_analysis/04_zscore_fromGSE.R)).


---
*Note: this README file has been generated automatically.* <br>
*Please do not modify it directly but instead work on [this config file](resources/config.yaml).*


