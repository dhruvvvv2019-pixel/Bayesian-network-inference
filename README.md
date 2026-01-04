# Bayesian-network-inference
Parameter learning and inference in Bayesian Networks for weather forecasting.

This project implements statistical learning techniques to estimate missing conditional probability tables (CPTs) in a large Bayesian Network modeling atmospheric conditions and severe weather outcomes. Using historical weather records with missing values, the system learns unknown parameters in the provided hailfinder.bif network and produces a complete Bayesian model suitable for probabilistic inference tasks such as hail and severe storm prediction.

The implementation adheres strictly to the constraints: no external Bayesian Network or EM libraries are used, and learning is completed efficiently on datasets with over 10,000 records.
