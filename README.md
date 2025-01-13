# SHP-1i Predictive Model Using Random Forest
**Method summary**

This R script implements the SHP-1i in-silico treatment predictive model using a Random Forest (RF) multi-class classification approach (see Figure 1). 
To develop this model, we first identified 121 differentially expressed genes (DEGs, FDR < 0.1) in human macrophages following in-vitro SHP-1i treatment. During model construction, these genes were then narrowed down to 119 through LASSO feature selection. Using the original, unaltered expression data for these 119 genes from 654 human plaque samples, we trained a random forest classification model to predict plaque subtypes. The data was partitioned into 1,000 randomly generated subsets, each split into 80% training and 20% testing sets. RF model performance was evaluated using multiclass area under the ROC curve (AUC), with an average AUC of 75% or higher considered indicative of a well performing model.  

To simulate in-silico treatment of SHP-1i, the expression changes of the same 119 DESeq2-derived, LASSO-selected DEGs were projected onto the 654 human plaque samples. This was achieved by multiplying the expression values by the antilogarithmic transformed fold changes, resulting in a “perturbed” dataset. The trained predictive model was then applied to classify the perturbed dataset, and and these classifications were used to quantify the effects of the in-silico SHP-1i treatment. After 1,000 iterations,  an averaged output was generated, revealing predicted changes in plaque classification after SHP-1i in-silico treatment.   

The shifts in plaque subtypes following this projection were assessed for statistical significance using McNemar’s test, which evaluates paired data to determine whether the observed shifts (e.g., from more vulnerable to less vulnerable plaques) are meaningful. To further the interpret results, odds ratios (ORs) were calculated by dividing the number of plaques predicted to shift favorably (e.g., transitioning from more vulnerable to less vulnerable subtypes) by the number of plaques predicted to shift unfavorably (e.g., transitioning from less vulnerable to more vulnerable subtypes). Corresponding confidence intervals were calculated using standard formulas and p-values were derived from McNemar's Chi-squared test. 

Libraries used:
progress 1.2.3,
pROC 1.18.5,
caret 6.0-94,
glmnet 4.1-8,
randomForest 4.7-1.1,
EnsDb.Hsapiens.v86 2.99.0,
biomaRT 2.60.1,
ggalluvial 0.12.5,
limma 3.60.4,
ggplot2 3.5.1
dplyr 1.1.4

Figure 1 Overview of model building 
![VAR flow chart-4](https://github.com/user-attachments/assets/0cf59108-5537-429f-95d7-3ff7deebbc78)
Numbered boxes indicate corresponding phase in step-by-step description. 

## Getting started: step by step description 
Before initiating the predictive Random Forest model loop, two main dataframes need to be prepared:

### 1. Before treatment data frame 
This data frame will be used to **train and test** the RF model. It contains the original normalized bulk RNA expression of the patient population, which will help the model learn patterns to predict outcomes. Additionally, it includes the known plaque subtype for each patient.

### 2. After in-silico treatment data frame
This data frame represents data after simulated (in-silico) treatment. It will be **fed into the trained Random Forest model** to make treatment predictions based on the learned patterns from the training phase. It **does not** contain any plaque subtype information.


## Steps to prepare data frames (script 1)

1. **Load and normalize data**:  
   The raw bulk RNA dataset is loaded, scaled, and quantile normalized using the `limma` package (with patients as rows and genes as columns).

2. **Load in-vitro cell drug treatment data**:  
   The in-vitro cell drug treatment expression data (for the drug of interest or controls) is loaded and saved in a table with the following variables: `gene` (ENSG), `symbol`, and `FoldChange`. Genes that are not differentially expressed after treatment (FRD>0.1) are filtered out. The variable `FoldChange` is created by transforming the existing log fold change values using the formula:
   
   $$
   \text{FoldChange} = 2^{\text{log fold change}}
   $$

3. **Filter normalized bulk RNA data frame**:  
   The normalized bulk RNA data frame is filtered for plaque subtype-defining genes (No. = 5,000) and genes that were differentially expressed after the in-vitro drug experiment (FDR < 0.1). Additionally, genes with low counts (count <= 10) are excluded.

4. **Create after in-silico treatment data frame**:  
   The in-silico treatment data frame is created by multiplying the expression values from the 'Before Treatment DataFrame' with the corresponding gene 'FoldChange' values from the in-vitro drug experiment table.

## RF model building and performance evaluation (script 2)

5. **Initialize RF Prediction model loop**:  
   The RF model is trained and tested using the 'Before Treatment DataFrame' (See Figure 2). Model performance is evaluated using the Area Under the Curve (AUC), and the predicted plaque subtype is calculated by feeding the trained RF model with the 'After In-Silico Treatment DataFrame.' To overcome overfitting, the RF prediction model is performed 1,000 times, with each iteration utilizing a different subset of the patient population created using `createDataPartition(y=..., p=0.80, times = 1000, list = TRUE)`. This results in a final output that represents the average prediction over 1,000 separate models.

Figure 2: Zoomed-in overview of model building and performance evaluation 
![Var model building and QC-5](https://github.com/user-attachments/assets/b393c3b8-ef53-4794-8484-aece84a58182)

## Model results & visualization (script 3)
6. Differences in the proportion of plaque subtype before and after in-silico treatment were assessed using Chi-sqaured test, with p-value<0.05 concidered signifciant. Differences between simplified grouping "more vulnerable" or "less vulnerable" plaques were determined using McNemar's chi-squared test, with p-value<0.05 concidered signficant. Model results were visualized in an alluvial plot using ggalluvial() R function.  

Formulas used to calculate Odds ratios and corresponding upper and lower CIs:
   a= No. less vulnerable plaques after in-silico treatment
   b= No. more vulnerable plaques  after in-silico treatment
   
   $$
   \text{Odds Ratio} = \frac{a}{b}
   $$
 
   $$
   \text{SE}_{\log(\text{OR})} = \sqrt{\frac{1}{a} + \frac{1}{b}}
   $$

   $$
   \log(\text{OR}) = \log(\text{Odds Ratio})
   $$
   
   $$
   \text{Lower CI} = \exp(\log(\text{OR}) - 1.96 \times \text{SE}_{\log(\text{OR})})
   $$
  
   $$
   \text{Upper CI} = \exp(\log(\text{OR}) + 1.96 \times \text{SE}_{\log(\text{OR})})
   $$

# RF Output

![image](https://github.com/user-attachments/assets/3faa77ec-978b-4137-af16-4f63c2fe438e)

![image](https://github.com/user-attachments/assets/c488dafd-de8f-44c5-a66e-aa2bb9656cf1)

![image](https://github.com/user-attachments/assets/ca400bd0-5556-4d4b-91dc-add4431b6359)







