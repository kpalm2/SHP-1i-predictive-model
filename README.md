# SHP-1i Predictive Model Using Random Forest

This R script implements the SHP-1i in-silico treatment predictive model using a Random Forest multi-class classification approach (see flow-chart).

To do this, we projected DESeq2 derived 119 differentially expressed genes (FDR<0.1), further selected via lasso feature selection, that occur in human macrophages following in-vitro SHP-1i onto data from 654 human plaque samples. The processed data was divided into 1000 randomly partitioned subsets, with each subset split into 80% training and 20% testing sets to predict plaque subtypes based on predicted transcriptomic changes and to overcome overfitting. The resulting trained models were then evaluated using multiclass area under the ROC curve (AUC), with an average AUC of 75% or higher considered indicative of a well performing model. From these models, an averaged output from 1000 iterations was generated, revealing predicted changes in plaque classification after SHP-1i in-silico treatment. Significant predicted shifts in plaque phenotypes after in-silico treatment were determined using Chi-sqaured test or the McNemar's chi-squared test with p<0.05 considered significant. Favorable shift odds ratios were manually calculated by dividing the number of least vulnerable plaques after in-silico treatment by the number of most vulnerable plaques after in-silico treatment.Corresponding confidence intervals were calculated using standard formulas (ref) and p-values were derived from Mcnemars Chi-squared test.

**Flow chart methods**
zoomed out flow-chart TBA
training flow chart TBA
Minimal gene seleciton graph TBA



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

## Getting Started

Before initiating the predictive Random Forest model loop, two main dataframes need to be prepared:

### 1. Before Treatment DataFrame
This dataframe will be used to **train and test** the Random Forest model. It contains the original normalized bulk RNA expression of the patient population, which will help the model learn patterns to predict outcomes. Additionally, it includes the known plaque subtype for each patient.

### 2. After In-Silico Treatment DataFrame
This dataframe represents data after simulated (in-silico) treatment. It will be **fed into the trained Random Forest model** to make treatment predictions based on the learned patterns from the training phase. It **does not** contain any plaque subtype information.

## Steps to Prepare Dataframes

1. **Load and Normalize Data**:  
   The raw bulk RNA dataset is loaded, scaled, and quantile normalized using the `limma` package (with patients as rows and genes as columns).

2. **Load In-Vitro Cell Drug Treatment Data**:  
   The in-vitro cell drug treatment expression data (for the drug of interest or controls) is loaded and saved in a table with the following variables: `gene` (ENSG), `symbol`, and `FoldChange`. Genes that are not differentially expressed after treatment (FRD>0.05) are filtered out. The variable `FoldChange` is created by transforming the existing log fold change values using the formula:

   $$
   \text{FoldChange} = 2^{\text{log fold change}}
   $$

3. **Filter Normalized Bulk RNA Dataframe**:  
   The normalized bulk RNA dataframe is filtered for plaque subtype-defining genes (n = 5000) and genes that were differentially expressed after the in-vitro drug experiment (FDR < 0.05). Additionally, genes with low counts (count >= 10) are excluded.

4. **Create After In-Silico Treatment DataFrame**:  
   The in-silico treatment dataframe is created by multiplying the expression values from the 'Before Treatment DataFrame' with the corresponding gene FoldChange values from the in-vitro drug experiment table.

5. **Initialize Random Forest Prediction Model Loop**:  
   The model is trained and tested using the 'Before Treatment DataFrame.' Model performance is evaluated using the Area Under the Curve (AUC), and the predicted plaque subtype is calculated by feeding the trained model with the 'After In-Silico Treatment DataFrame.' The Random Forest prediction model is performed 1000 times, with each iteration utilizing a different subset of the patient population created using `createDataPartition(y=..., p=0.80, times = 1000, list = TRUE)`. This results in a final output that represents the average prediction over 1000 separate models.

6. **Statistical analysis of results**:
   Differences in the proportion of plaque subtype before and after in-silico treatment were assessed using Chi-sqaured test, with p-value<0.05 concidered signifciant. Differences between simplified grouping "more vulnerable" or "less vulnerable" plaques were determined using McNemar's chi-squared test, with p-value<0.05 concidered signficant.

  **Formulas used to calculate Odds ratios and corresponding upper and lower CIs**

   Odds ratio
   a= No. less vulnerable plaques after in-silico treatment
   b= No. more vulnerable plaques  after in-silico treatment
   
   $$
   \text{Odds Ratio} = \frac{a}{b}
   $$

   Standard error and 95% CI for log odds ratio
   
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

# Example Output

![image](https://github.com/user-attachments/assets/3faa77ec-978b-4137-af16-4f63c2fe438e)

![image](https://github.com/user-attachments/assets/c488dafd-de8f-44c5-a66e-aa2bb9656cf1)

![image](https://github.com/user-attachments/assets/ca400bd0-5556-4d4b-91dc-add4431b6359)





