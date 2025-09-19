import os

cuts = [20, 35, 55]
hws = [5]
dataset_names = ["Sgr A", "Sgr B2", "Sgr C"]
actual_names = ["SGRA", "SGRB", "SGRC"]

for i, dataset in enumerate(actual_names):
    for cut in cuts:
        for hw in hws:
            #low_command = "python /rascil-main/rascil/apps/rascil_imager.py --ingest_msname " + dataset + "_small_baselines_" + str(cut + hw) + \
            #    ".ms --ingest_vis_nchan 1 --imaging_npixel 512 --imaging_cellsize 0.00001849451 --imaging_weighting uniform " + \
            #    "--clean_nmajor 5 --clean_algorithm mstep --clean_fractional_threshold 0.3 --clean_threshold 1e-3 --clean_restored_output integrated " + \
            #    "--mstep_output_intermediate True --mstep_mode low --mstep_lambda 0.05 --mstep_lambda_mul 2 --mstep_wavelet daubechies " + \
            #    "--mstep_cut_center " + str(cut) + " --mstep_cut_hw " + str(hw)

            low_path = "/sunwang/radio-imaging/results/cut_results/" + dataset + "/" + str(cut) + "/" + str(hw) + "hw/low/"
            multistep_path = "/sunwang/radio-imaging/results/cut_results/" + dataset + "/" + str(cut) + "/" + str(hw) + "hw/high/"

            multistep_command = "python /rascil-main/rascil/apps/rascil_imager.py --ingest_msname " + dataset + "_long_baselines_" + str(cut - hw) + \
            ".ms --ingest_vis_nchan 1 --imaging_npixel 512 --imaging_cellsize 0.00001849451 --imaging_weighting uniform " + \
            "--clean_nmajor 5 --clean_algorithm mstep --clean_fractional_threshold 0.3 --clean_threshold 1e-3 --clean_restored_output integrated " + \
            "--mstep_output_intermediate True --mstep_mode multi-step --mstep_lambda 0.05 " + \
            "--mstep_lambda_mul 2 --mstep_lowfreq_image " + low_path + dataset + "_small_baselines_" + str(cut + hw) + "_nmoment4_cip_deconvolved.fits " + \
            " --mstep_wavelet daubechies --mstep_cut_center " + str(cut) + " --mstep_cut_hw " + str(hw)

            
            #print("RUNNING: " + low_command)
            #os.system(low_command)
            #os.system("rm " + low_path + "*")
            #os.system("cp *.fits " + low_path)
            print("RUNNING: " + multistep_command)
            os.system(multistep_command)
            os.system("rm " + multistep_path + "*")
            os.system("cp *.fits " + multistep_path)
            os.system("rm *.fits")