# Calculate average entropy ratio for denoising method (average across all sites).
import os, re
from statistics import mean

def calc_ent_ratios(data_dir:str, ent_options:list=["ent", "no_ent"], alpha_options:list=[1, 3, 5, 7, 9, 11, 13]):
    """
    
    """

    avg_ent_ratios = {}

    for y_n in ent_options:
        for alpha in alpha_options:
            path_to_logs = f"{data_dir}/benchmarking/{y_n}_alpha_{alpha}_denoised/logs"
            ent_ratios = {}

            for log in os.listdir(path_to_logs):

                with open(f"{path_to_logs}/{log}") as file:
                    lines = file.readlines()
                    for line in lines:
                        if "second" in line:
                            num_match = re.findall(pattern="\d+\.\d+", string=line)
                            if num_match:
                                ent_2 = float(num_match[0])
                                print(ent_2)
                        if "third" in line:
                            num_match = re.findall(pattern="\d+\.\d+", string=line)
                            if num_match:
                                ent_3 = float(num_match[0])
                                print(ent_3)
                        
                    ent_ratios[log] =  ent_2 / ent_3
            #avg_ent_ratios[f"{y_n}_{alpha}"] = mean(list(ent_ratios.values()))

    print(ent_ratios)


if __name__ == "__main__":
    calc_ent_ratios(data_dir="../../data/test_data")