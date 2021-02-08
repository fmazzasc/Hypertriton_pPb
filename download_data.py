import os

if not os.path.exists("Trees"):
        os.makedirs("Trees")

if not os.path.exists("Tables"):
        os.makedirs("Tables")

if not os.path.exists("Utils"):
        os.makedirs("Utils")
if not os.path.exists("Utils/ProdModels"):
        os.makedirs("Utils/ProdModels")

if not os.path.exists("Trees/HyperTritonTree_13bc.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/2Body/HyperTritonTree_13bc.root  Trees/.")

if not os.path.exists("Trees/HyperTritonTree_13bc_LS.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/2Body/HyperTritonTree_13bc_LS.root  Trees/.")

if not os.path.exists("Trees/HyperTritonTree_16qt.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/2Body/HyperTritonTree_16qt.root  Trees/.")

if not os.path.exists("Trees/HyperTritonTree_16qt_LS.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/2Body/HyperTritonTree_16qt_LS.root  Trees/.")

if not os.path.exists("Trees/HyperTritonTree_17d.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/2Body/HyperTritonTree_17d.root  Trees/.")    ##G3 based MC

if not os.path.exists("Trees/HyperTritonTree_20l2.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/2Body/HyperTritonTree_20l2.root  Trees/.")   #G4 based MC

if not os.path.exists("Utils/AnalysisResults_pPb.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/2Body/AnalysisResults_pPb.root  Utils/.")

if not os.path.exists("Utils/ProdModels/s3_3body.csv"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/pPbProdModels/s3_3body.csv  Utils/ProdModels/.")

if not os.path.exists("Utils/ProdModels/s3_2body.csv"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/pPbProdModels/s3_2body.csv  Utils/ProdModels/.")

if not os.path.exists("Utils/ProdModels/FullCSM-S3.dat"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/pPbProdModels/FullCSM-S3.dat  Utils/ProdModels/.")

if not os.path.exists("Utils/ProdModels/CSM_predictions_S3_T155MeV_Vc3dNdy.root"):
    os.system("scp lxplus.cern.ch:/eos/user/h/hypertriton/trees/pPbProdModels/CSM_predictions_S3_T155MeV_Vc3dNdy.root  Utils/ProdModels/.")

if not os.path.exists("Utils/AbsorptionHe3/"):
    os.system("scp -r lxplus.cern.ch:/eos/user/h/hypertriton/trees/AbsorptionHe3/  Utils/.")

