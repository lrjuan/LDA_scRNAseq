## LDA_scRNAseq

Source codes and test data for the LDA-based scRNA-seq data analysis workflow.

The workflow consists of a series tools, primarily aimed at applying the LDA model to single-cell RNA sequencing datasets, and contructing a cell-function-gene three-layer framework to help researchers better understanding data.

### The LDA model

LDA is the abbreviation of Latent Dirichlet Allocation. LDA is a probabilistic topic model utilizing unsupervised learning that was initially proposed in Natural Language Processing (NLP) area. The basic assumption is that the reason we observe a specific set of words in a document is actually determined by a group of latent attributes in the document, i.e., topics.

The LDA model was originally developed for mining text from large-scale corpora, aiming for latent ‘topic’ recognition in massive observed documents. In the LDA model, each document is regarded as a mixture of multiple topics, each topic associates to lots of words. Based on the above assumption, specific words were chosen to be written in a document, is because the document focuses on the associated topics.

More generally, the LDA model is applied to identify latent attributes that are difficult to directly observe from data. These datasets usually have two-layer structures, including data units and their unordered collections. The observed collection-unit relationship pattern is largely determined by the latent attributes. In many cases, recognizing the latent attributes is the main purpose of such a data analysis process.

In the scRNA-seq analysis context, **cell** is regarded as '**document**', **gene** is regarded as '**word**', ***putative function (PF)*** is regarded as '**topic**'.

## Usage

### List of tools and example data

* run_lda.pl - main program
* cell_cluster_annotation.R - downstream analysis tools based on the LDA results
* main.R - example commands of downstream analysis tools
* test.txt(.gz)  - example data, simulated by Splatter
* test.labels.txt - the group labels of the test data
* run_lda.mex.pl - main program for CellRanger Count format: MEX
* transpose.pl - transpose a gene x cell matrix to cell x gene matrix

### run_lda.pl

run_lda.pl is the main program of this workflow. It takes the cell-gene expression matrix as input, extracts PFs (topics) from the data, generates PF distributions for each cell and Gene distributions for each PF. 

run_lda.pl is an all-in-one style PERL script, the actual LDA modeling is accomplished by MALLET. To use run_lda.pl, **MALLET is required**.

>**[MALLET](https://mimno.github.io/Mallet/)** stands for **MA**chine **L**earning for **L**anguag**E** **T**oolkit, is a Java-based package for statistical natural language processing, document classification, clustering, topic modeling, information extraction, and other machine learning applications to text.

The MALLET package and installation instruction can be found at [https://github.com/mimno/Mallet](https://github.com/mimno/Mallet).

Assuming the MALLET has been deployed in the system. Before running run_lda.pl, several global parameters neet to be configured.

|Global parameter|Value type|Default|Note|
|:----|:----|:----|:----|
|Threads_num|INTEGER|1|How many threads can be occupied during the LDA training.|
|Mallet_path|DIRECTORY|/Mallet/bin|The pre-installed Mallet program location.|
|Iterations_num|INTEGER|500|Number of iterations for Gibbs sampling.|

The global parameters can be found at the top of the PERL script.

Once the global parameters were configured properly accroding to the current system, everything is ready.

Usage:  `perl run_lda.pl -input Path_of_expression_matrix -output Prefix_of_results -topics Number_of_topics`

The program recieves 3 parameters: -input, -output, and -topics

|Parameter|Value type|Note|
|:----|:----|:----|
|-input|FILE|Single input file: **Cells as rows, genes as columns**. The first row is gene ID, the first column is cell ID.|
|-output|FILENAME|The results include 5 files:<br>(1) prefix.cell.txt - Cell list<br>(2) prefix.gene.txt - Gene list<br>(3) prefix.c2k.txt - Cell-PF(Document-Topic) matrix<br>(4) prefix.k2g.txt - PF-Gene(Topic-Term) matrix<br>(5) prefix.mallet.log - Log for the Mallet run|
|-topics|INTEGER|Topics number for LDA training|

All parameters are required.

In the command line, users can use the --help/-h parameter to view the usage of the program.

`perl run_lda.pl --help`

### cell_cluster_annotation.R

cell_cluster_annotation.R is the cluster annotation and function interpretation methods. It contains 6 functions, including read_cellTopic, read_topicGene, LDAFindClusters, clusterHGT, representTopicCluster, and representTopicCluster. It takes the results of the LDA modelling and produces functional interpretations, such as the reading of LDA result, cell clusters and their annotations, representative genes and putative functions, etc.

#### read_cellTopic

read_cellTopic is a function to read Cell-PF(Document-Topic) matrix.

Usage:  `read_cellTopic(outPrefix)`

| Parameter | Value type | Note                           |
| :-------- | :--------- | :----------------------------- |
| outPrefix | FILENAME   | Path and Prefix of LDA results |

#### read_topicGene

read_topicGene is a function to read PF-Gene(Topic-Term) matrix.

Usage:  `read_topicGene(outPrefix)`

| Parameter | Value type | Note                           |
| :-------- | :--------- | :----------------------------- |
| outPrefix | FILENAME   | Path and Prefix of LDA results |

#### LDAFindClusters

LDAFindClusters is a function to identify clusters of cells by LDA and Hellinger distance based K-medoids algorithm

Usage:  `LDAFindClusters(cellTopic, clusterNumber)`

| Parameter     | Value type | Note                           |
| :------------ | :--------- | :----------------------------- |
| cellTopic     | MATRIX     | Cell-PF(Document-Topic) matrix |
| clusterNumber | INTEGER    | specify the number of clusters |

#### clusterHGT

clusterHGT is a function to automatically predict cell type by identifying representative genes for each cell clusters and performing hypergeometric test against pre-established cell type marker LISTs.

Usage:  `clusterHGT(cellTopic, topicGene, clusterResult, pathways , n.features , minSize = 0, log.trans = TRUE, p.adjust = TRUE)`

| Parameter     | Value type | Default | Note                                                         |
| :------------ | :--------- | :------ | :----------------------------------------------------------- |
| cellTopic     | MATRIX     |         | Cell-PF(Document-Topic) matrix                               |
| topicGene     | MATRIX     |         | PF-Gene(Topic-Term) matrix                                   |
| clusterResult | VECTOR     |         | clustering result                                            |
| pathways      | LIST       |         | pre-established cell type marker lists                       |
| n.features    | INTEGER    |         | top n representative genes to consider for hypergeometric test |
| minSize       | INTEGER    | 0       | minimum number of overlapping genes in pathways              |
| log.trans     | BOOLEAN    | TRUE    | if TRUE tranform the pvalue matrix with -log10 and convert it to sparse matrix |
| p.adjust      | BOOLEAN    | TRUE    | if TRUE apply Benjamini Hochberg correction to p-value       |

#### representTopicCluster

representTopicCluster is a function to identify representative PFs for cell clusters.

Usage:  ` representTopicCluster(cellTopic, clusterResult)`

| Parameter     | Value type | Note                           |
| :------------ | :--------- | :----------------------------- |
| cellTopic     | MATRIX     | Cell-PF(Document-Topic) matrix |
| clusterResult | VECTOR     | clustering result              |

#### representTopicCluster

representTopicCluster is a function to identify representative genes for cell clusters.

Usage:  ` representGeneCluster(cellTopic, topicGene, clusterResult, n.features)`

| Parameter     | Value type | Note                                                         |
| :------------ | :--------- | :----------------------------------------------------------- |
| cellTopic     | MATRIX     | Cell-PF(Document-Topic) matrix                               |
| topicGene     | MATRIX     | PF-Gene(Topic-Term) matrix                                   |
| clusterNumber | INTEGER    | specify the number of clusters                               |
| n.features    | INTEGER    | specify how many top representative genes will be extracted. |

### main.R

main.R is an example commands of cell_cluster_annotation.R. It shows how to process the test data using cell_cluster_annotation.R. Its input is LDA result(PF = 25). This input is generated by processing cell-gene matrix of test data with run_lda.pl. The command to generate this input is ` perl run_lda.pl -input testData/cell_gene.txt -output testData/PF25 -topics 25`

### test.tar.gz

The test data is sampled from the PBMC dataset (Zheng, et al., 2017), with 500 samples taken from each cell type, for a total of 5,000 cells. Each cell has 21,758 genes.

There are 8 files in the compressed packages. 

|File|Description|
|:----|:----|
|cell_gene.txt|cell x gene expression matrix|
|cell_labels.csv|the corresponding cell labels|
|cell_types_marker.csv|gene markers of the cell types|
|PF25.c2k.txt|LDA modelling results (PF = 25)|
|PF25.cell.txt|LDA modelling results (PF = 25)|
|PF25.gene.txt|LDA modelling results (PF = 25)|
|PF25.mallet.log|LDA modelling results (PF = 25)|

### test.txt(.gz) and test.labels.txt

**test.txt** is an example data for run_lda.pl. The file is a cell x gene matrix, cells as rows, genes as columns. The first row is gene ID, the first column is cell ID.

There are 1,000 cells in the test data, each has 15,000 genes. The expression data is simulated by Splatter. 

The 1,000 cells can be grouped into 5 groups. **test.labels.txt** records the correspondence between cell IDs and their grouping.

### run_lda.mex.pl

run_lda.mex.pl is an alternative version of main program. Everything is same to run_lda.pl except it takes directory of the 'CellRanger Count' results as input parameter.

Usage:  `perl run_lda.mex.pl -input Path_of_MEX_directory -output Prefix_of_results -topics Number_of_topics`

|Parameter|Value type|Note|
|:----|:----|:----|
|-input|DIRECTORY|The Path_of_MEX_directory is the "CellRanger Count" output direcotry, must contain 3 files:<br>(1) features.tsv(.gz) or genes.tsv(.gz)<br>(2) barcodes.tsv(.gz) or cells.tsv(.gz)<br>(3) matrix.mtx(.gz)<br>*The above files can be kept gzipped.|
|-output|FILENAME|The results include 3 files:<br>(1) prefix.c2k.txt - Cell-PF(Document-Topic) matrix<br>(2) prefix.k2g.txt - PF-Gene(Topic-Term) matrix<br>(3) prefix.mallet.log - Log for the Mallet run<br>*run_lda.mex.pl does not provide cell list and gene list, because they exist independently in the input directory in the first place.|
|-topics|INTEGER|Topics number for LDA training|

An example data for run_lda.mex.pl can be found [Here](https://cf.10xgenomics.com/samples/cell-exp/7.0.0/1k_mouse_kidney_CNIK_3pv3/1k_mouse_kidney_CNIK_3pv3_filtered_feature_bc_matrix.tar.gz)

### transpose.pl

run_lda.pl takes **cell x gene** matrix as input, i.e., cells as rows, genes as columns. This is because files are stored line by line on the hard disk. **Cell x gene** matrix enables the program to load the entire gene expression of a cell each time, rather than loading the expression of a gene in all cells each time. 

However, many softwares, including R and R packages, store and consider matrix objects as a column-by-column vector. Thus a lot of scRNA-seq tools take **gene x cell** matrix as input, or produce such matrix as output.

Here we provide transpose.pl program to transpose the gene x cell matrix into cell x gene matrix, to meet the input format requirement of run_lda.pl.

Usage:  `perl transpose.pl Path_of_input_file > Output_file`

Memory consumption for converting large-scale data is high, so please be cautious while running it.

---

Please do not hesitate to address comments/questions/suggestions regarding this tool to: pgbrowser@gmail.com.

L.J. Apr. 8, 2023
