# DPREx
DPREx is a machine-learning based tool to discriminate potential genes intervals to be pathogenic in which repeats
mutation happens.

## how does the model work?
First, we use comprehensive epigenetic data which include chromatin accessibility, histone modification marks, the
binding strength of TFs, distance to alternative splicing, and non-B DNA structure to annotate the intervals of input
repeats.

![workflow](./workflow.png)

Then, a pre-trained XGBoost based model is loaded to predict the annotated table. This model uses features dominantly
come from repeat interval annotation, and doesn't care about the detail length of the interval , making its prediction
mainly based on gene expression and regulation, instead of pattern of mutations.

## how to use
* prepare for the environment
    * make sure the python's version is not less than 3.8. 
    * use virtualenv to create a new environment (optional but recommended).
    * type in: `install -r requirement.txt`, make sure you are in the same path of this project.

* config epigenetic data resources in manifest.txt.
    * It should contain your sample file and epigenetic data file (in various formats can be found in the sample.).
    * If you want to add new property, add new line following.
    * To remove property, use '%' to comment out that line.

* set configuration in config.json.
* organise your input file in BED3 format.
* invoke the command line API:<br>
  `python main.py path-to-input/input.bed -o save_dir/`

## parameter
* input [positional] : the path of input BED6 file
* --output [-o] : the path to store the result csv
* --help <-h> : for help

## input file
The input file should be organised following the BED6 format, the columns are:<br>
chromosome start end<br>
seperated with tab

Here are some examples:<br>
chr1 12252260 12252264<br>
chr3 187388896 187388910<br>
chr4 1019055 1019075

## result file
* The result will be saved in a CSV file, with columns called repeats_id, probability_score and pathogenic. 
* For a record like "chr1 12252260 12252264", its repeats_id is "chr1@12252260@12252264". 
* The order of the records can be different from the input file, but it can be re-sorted with repeats_id.
* The pathogenic column denotes whether a mutation is pathogenic if it happens in the repeat interval, the
probability_score is positively related to that likelihood.
