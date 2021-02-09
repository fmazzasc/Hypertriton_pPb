# Hypertriton_pPb

Contains the code and the data to study the hypertriton production in pPb collisions at 5.02TeV with ALICE. Analysis based on xgboost Boosted Decision Trees Classifier.
## Run the analysis
- Download the data: `python download_data.py`
- Generate flat trees for feeding the BDT: `cd GenerateTables/`, then `python table_generator.py`
- Train and test the BDT model: `cd Analysis/`, then `python training_and_testing`. The file `Config.yaml` allows you to customize the training and define a BDT      efficiency range for computing the systematic variations
- Extract the hypertriton signal, compute the yield and the S3 value: `python compute_s3.py`
- Extract the hypertriton signal, compute the yield and the hyp/Lambda ratio: `python compute_hyp_lambda_ratio.py`
- Check the results: `cd Results/`
