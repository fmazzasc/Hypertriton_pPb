# Hypertriton_pPb

Contains the code and the data to study the hypertriton production in pPb collisions at 5.02TeV with ALICE. Analysis based on xgboost Boosted Decision Trees Classifier.
## Run the analysis
- download the data: `python download_data.py`
- generate flat trees for feeding the BDT: `cd GenerateTables/`, then `python table_generator.py`
- train and test the BDT model: `cd Analysis/`, then `python training_and_testing`. The file `Config.sh` allows you to customize the training and define a BDT      efficiency range for computing the systematic variations
- extract the hypertriton signal, compute the yield and the s3 value: `python compute_yield.py`
- Check the results: `cd Results/`
