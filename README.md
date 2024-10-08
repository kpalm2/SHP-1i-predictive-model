# SHP-1i Predictive Model Using Random Forest

This R script implements the SHP-1i in-silico treatment predictive model using a Random Forest multi-class classification approach.

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
   The in-vitro cell drug treatment expression data (for the drug of interest or controls) is loaded and saved in a table with the following variables: `gene` (ENSG), `symbol`, and `FoldChange`. genes that are not differentially expressed after treatment (FRD>0.05) are filtered out. The variable `FoldChange` is created by transforming the existing log fold change values using the formula:

   $$
   \text{FoldChange} = 2^{\text{log fold change}}
   $$

3. **Filter Normalized Bulk RNA Dataframe**:  
   The normalized bulk RNA dataframe is filtered for cluster-defining genes (n = 5000) and genes that were differentially expressed after the in-vitro drug experiment (FDR < 0.05). Additionally, genes with low counts (count >= 10) are excluded.

4. **Create After In-Silico Treatment DataFrame**:  
   The in-silico treatment dataframe is created by multiplying the expression values from the 'Before Treatment DataFrame' with the corresponding gene FoldChange values from the in-vitro drug experiment table.

5. **Initialize Random Forest Prediction Model Loop**:  
   In this loop, the model is trained and tested using the 'Before Treatment DataFrame.' Model performance is evaluated using the Area Under the Curve (AUC), and the predicted plaque subtype is calculated by feeding the trained model with the 'After In-Silico Treatment DataFrame.' The Random Forest prediction model is performed 1000 times, with each iteration utilizing a different subset of the patient population created using `createDataPartition(y=..., p=0.80, times = 1000, list = TRUE)`. This results in a final output that represents the average prediction over 1000 separate models.


# Example Output

![image](https://github.com/user-attachments/assets/3faa77ec-978b-4137-af16-4f63c2fe438e)

![image](https://github.com/user-attachments/assets/c488dafd-de8f-44c5-a66e-aa2bb9656cf1)

![image](https://github.com/user-attachments/assets/ca400bd0-5556-4d4b-91dc-add4431b6359)





