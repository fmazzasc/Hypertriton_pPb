import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
import hipe4ml.analysis_utils as au
import hipe4ml.plot_utils as pu
import mplhep as mpl
import xgboost as xgb
import helpers as hp
import yaml

matplotlib.use('pdf')
plt.style.use(mpl.style.ALICE)


with open("Config.yaml", 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

training = params['TRAINING']
bkg_fraction = params["BACKGROUND_OVER_SIGNAL"]
optmize = params['OPTIMIZE']
application = params['APPLICATION']
significance_scan = params["SIGNIFICANCE_SCAN"]
working_point = params['EFF_WORKING_POINT']
variation_range = params['SYST_VARIATION_RANGE']
data_table_name = params['DATA_TABLE']
bkg_table_name = params['BKG_TABLE']
signal_table_name = params['SIGNAL_TABLE']
MODEL_PARAMS = params['XGBOOST_PARAMS']
HYPERPARAMS = params['HYPERPARAMS']

path_to_data = "../Tables/"
results_ml_path = "../Results/PlotsML"
ml_model_path = "../Utils/Models"
utils_path = "../Utils"
efficiencies_path = "../Utils/Efficiencies"
selected_df_path = "../Utils/ReducedDataFrames"


print("---------------------------------------------")
print("Data loading...")

if training:

        signalH = TreeHandler(path_to_data + signal_table_name, "SignalTable")
        bkgH = TreeHandler(path_to_data + bkg_table_name, "DataTable")
        bkgH.get_data_frame().drop_duplicates(inplace=True)        

        if bkg_fraction!=None:
                bkgH.shuffle_data_frame(size=bkg_fraction*len(signalH), inplace=True, random_state=52)

        train_test_data = au.train_test_generator([signalH, bkgH], [1,0], test_size=0.5, random_state=42)


        training_columns = ['TPCnSigmaHe3','ct','V0CosPA','ProngsDCA','He3ProngPvDCA','PiProngPvDCA','He3ProngPvDCAXY','PiProngPvDCAXY','NpidClustersHe3','TPCnSigmaPi']

        if not os.path.exists(results_ml_path):
                os.makedirs(results_ml_path)

        distr = pu.plot_distr([bkgH, signalH], training_columns, bins=63, labels=['Signal',"Background"],colors=["blue","red"], log=True, density=True, figsize=(18, 13), alpha=0.3, grid=False)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(results_ml_path + "/features_distributions.png", bbox_inches='tight')
        corr = pu.plot_corr([signalH,bkgH], training_columns + ["m"], ['Signal',"Background"])
        corr[0].savefig(results_ml_path + "/correlations.png",bbox_inches='tight')

        print("---------------------------------------------")
        print("Data loaded. Training and testing ....")

        params_range = {
        "max_depth": (8, 18),
        "learning_rate": (0.07,0.15),
        "n_estimators": (150, 250),
        "gamma": (0.3,0.5),
        "min_child_weight": (3,8),
        "subsample": (0.5,1),
        "colsample_bytree": (0.3,1),
        }

        model_hdl = ModelHandler(xgb.XGBClassifier(), training_columns)
        model_hdl.set_model_params(MODEL_PARAMS)
        model_hdl.set_model_params(HYPERPARAMS)
        if optmize:
                model_hdl.optimize_params_bayes(train_test_data,params_range,'roc_auc',njobs=-1, init_points=10, n_iter=20)

        y_pred_test = model_hdl.train_test_model(train_test_data, True, True)

        bdt_out_plot = pu.plot_output_train_test(model_hdl, train_test_data, 100, True, ["Signal", "Background"], True, density=True)
        bdt_out_plot.savefig(results_ml_path + "/bdt_output.png")

        if not os.path.exists(ml_model_path):
                os.makedirs(ml_model_path)
        model_hdl.dump_model_handler(ml_model_path + "/model_hndl.pkl")

        eff_arr = np.round(np.arange(0.5,0.99,0.01),2)
        score_eff_arr = au.score_from_efficiency_array(train_test_data[3], y_pred_test, eff_arr)

        if not os.path.exists(efficiencies_path):
                os.makedirs(efficiencies_path)
        np.save(efficiencies_path + "/efficiency_arr.npy", eff_arr)
        np.save(efficiencies_path + "/score_efficiency_arr.npy",score_eff_arr)

        print("---------------------------------------------")
        print("Training done")

if application:

        print("---------------------------------------------")
        print("Starting application: ..")
        dataH = TreeHandler(path_to_data + data_table_name, "DataTable")
        signalH = TreeHandler(path_to_data + signal_table_name, "SignalTable")
        lsH = TreeHandler(path_to_data + bkg_table_name, "DataTable")
        lsH.get_data_frame().drop_duplicates(inplace=True)

        simH = TreeHandler(path_to_data + signal_table_name, "GenTable")

        presel_eff = len(signalH)/len(simH)
        bdt_eff_arr = np.load(efficiencies_path + "/efficiency_arr.npy")
        score_eff_arr = np.load(efficiencies_path + "/score_efficiency_arr.npy")

        model_hdl = ModelHandler()
        model_hdl.load_model_handler(ml_model_path + "/model_hndl.pkl")

        dataH.apply_model_handler(model_hdl)
        lsH.apply_model_handler(model_hdl)

        sign_plot = hp.significance_scan(bdt_eff_arr, score_eff_arr, dataH, presel_eff, working_point, variation_range)

        if significance_scan:
                plt.ylim((-0.3,5))
                sign_plot.savefig(results_ml_path + "/significance_scan.png")
                
        syst_mask = np.logical_and(bdt_eff_arr >= working_point - variation_range, bdt_eff_arr <= working_point + variation_range)
        bdt_eff_syst_arr = bdt_eff_arr[syst_mask]
        score_eff_syst_arr = score_eff_arr[syst_mask]

        selected_dataH = dataH.get_subset(f"model_output>{score_eff_syst_arr[-1]}")
        selected_lsH = lsH.get_subset(f"model_output>{score_eff_syst_arr[-1]}")


        np.save(efficiencies_path + "/bdt_eff_syst_arr.npy", bdt_eff_syst_arr)
        np.save(efficiencies_path + "/score_eff_syst_arr.npy", score_eff_syst_arr)

        if not os.path.exists(selected_df_path):
                os.makedirs(selected_df_path)
        selected_dataH.write_df_to_parquet_files(selected_df_path + "/selected_df_data")
        selected_lsH.write_df_to_parquet_files(selected_df_path + "/selected_df_ls")
  
        print("---------------------------------------------")
        print("Application done.")




