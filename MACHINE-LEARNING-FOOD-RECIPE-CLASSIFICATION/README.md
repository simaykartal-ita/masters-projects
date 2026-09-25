# Food Recipe Classification

This project was completed as part of a **Machine Learning course** and focuses on a binary classification problem involving recipe provenance.

The objective was to classify recipes as either **American** or **Italian** based on ingredient-related features represented by **TF-IDF scores**.

The dataset contained **4,934 recipes**, including:

- **3,200 labeled recipes** for model training
- **1,734 unlabeled recipes** for final prediction
- **40 ingredient-based features**

The target classes were:

- `1` → American
- `2` → Italian

## Methodology

Several supervised machine learning models were evaluated using **stratified 10-fold cross-validation**.

The final solution combined three models:

- RBF Support Vector Machine
- Random Forest
- Extra Trees

These models were combined using a **soft-voting ensemble**, allowing the final prediction to incorporate probability estimates from all three classifiers.

The SVM model included feature scaling through a pipeline, while the tree-based models were trained directly on the TF-IDF features.

The main evaluation metric was **classification accuracy**, which was converted into an estimated number of misclassified observations on the 1,734-recipe test set.

## Final Model

The final model used:

```text
Soft Voting Ensemble
├── RBF-SVM
├── Random Forest
└── Extra Trees
```

After cross-validation, the ensemble was trained on the full labeled dataset and used to predict the classes of the 1,734 unlabeled recipes.

The final predictions were exported as:

```text
y_pred.txt
```

with one predicted class (`1` or `2`) per line.

## Technologies

- Python
- Pandas
- NumPy
- Scikit-learn
- Support Vector Machines
- Random Forest
- Extra Trees
- Ensemble Learning
- Cross-Validation

## Notes

The original training and test datasets are not included in this repository.

This project demonstrates the use of model validation, ensemble learning, and reproducible machine learning workflows for a structured binary classification problem.
