# Code for the paper: **Machine learning-based DNA methylation score for fetal exposure to maternal smoking: development and validation in samples collected from adolescents and adults**
Authors: Sebastian Rauschert, Phillip E Melton, Anni Heiskala, Ville Karhunen, Graham Burdge, Jeffrey M Craig, Keith M Godfrey, Karen Lillycrop, Trevor A Mori, Lawrence J Beilin, Wendy H Oddy, Craig Pennell, Marjo-Riitta Järvelin, Sylvain Sebert, Rae-Chi Huang
*(in Press DOI 10.1289/EHP6076)*

This repository contains the code that was used to train the following four models for the development of an in utero smoke exposure score:

1. Random Forest
2. Gradient Boosting Machine
3. Support Vector Machine
4. Elastic Net Regression

The folder `code` contains the training and testing code, whereas the folder `models` contains the final gradient boosting machine, random forest and support vector machine `caret` objects to predict on
new, unseen data.

Further, the `code` folder contains code to create the Reese et al. score and the Richmond et al. score, which were the "gold standard" models we tested against in our
study.

Further instructions on how to use this repository can be found in the manuscript and in **Supplement 2** of the manuscript