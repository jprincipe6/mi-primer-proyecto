#!/bin/bash python3

#import pyplot
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

class Record():
    def __init__(self):
        self.stages = ["stage1", "stage2", "stage3", "stage4"]
        self.current_stage = ""
        self.stats_normal = {"stage1": [], "stage2": [], "stage3": [], "stage4"
:[]        }
        self.stats_custom = {"stage1": [], "stage2": [], "stage3": [], "stage4"
:[]        }
   
        
    def check_stage(self, line):
        if ("STAGE - 1") in line :
            self.current_stage = "stage1"
        if ("STAGE - 2") in line :
            self.current_stage = "stage2"
        if ("STAGE - 3") in line :
            self.current_stage = "stage3"
        if ("STAGE - 4") in line :
            self.current_stage = "stage4"

    def record_stat(self, line):
        if ("Normal" in line):
            tokens = line.rstrip().split(" ")
            data = tokens[-1]
            if ("%" in data):
                data = data[:-1]
            self.stats_normal[self.current_stage].append(data)

        if ("Custom" in line):
            tokens = line.rstrip().split(" ")
            data = tokens[-1]
            if ("%" in data):
                data = data[:-1]
            self.stats_custom[self.current_stage].append(data)

    def print_stats(self):
        for s in self.stages:
            print(s)
            print(self.stats_normal[s] )
            print(self.stats_custom[s] )       


    def save_bar_plots(self, positions=[0,1]):
        barwidth = 0.3
	data = positions
	#data = [0,len(self.stats_normal["stage1"]) -1 ]
        for i in data:
            ys_normal = []
            ys_custom = []
            for s in self.stages:
                if (s == "stage4"):
                    ys_normal.extend(self.stats_normal[s][i*3:i*3+3])
                    ys_custom.extend(self.stats_custom[s][i*3:i*3+3])
                else:
                    ys_normal.append(self.stats_normal[s][i])
                    ys_custom.append(self.stats_custom[s][i])
                    
            ys_normal = np.array(ys_normal).astype(int)
            ys_custom = np.array(ys_custom).astype(int)
            len_ys = len(ys_normal)
            print("ys_custom: ", ys_custom)
            xs_normal = range(1, len_ys)
            r1 = np.arange(len(ys_normal))
            r2 = [x + barwidth for x in r1]
            plt.xticks([r + barwidth for r in range(len(ys_normal))] , ["S1", "S2","S3","S4.1","S4.2","S4.3"])
            plt.bar(r1, ys_normal, width=barwidth, color="b", label="normal")
            plt.bar(r2, ys_custom, width=barwidth, color="r", label="custom")
            plt.legend()

            # value of bar as label for each bar
            for j, v in enumerate(ys_custom):
                    plt.text(r2[j] - 0.05, v + 0.25, str(v))
 
            for j, v in enumerate(ys_normal):
                    plt.text(r1[j] - 0.1, v + 0.25, str(v))
                                
            plt.grid(axis="y", linestyle="-")
            plt.savefig("stats_cycle_"+str(i)+".png")
            plt.close()

    def save_data_plots(self):
        width = 0.3
        for s in self.stages:
            len_ys= len(self.stats_normal[s])
            if (len_ys > 0):
                ys_normal = self.stats_normal[s]
                ys_custom = self.stats_custom[s]
                xs_normal = range(0, len_ys)
                xs_custom = np.arange(0, len_ys, 1)
                fig, ax = plt.subplots()
                plt.gca().invert_yaxis()
                ax.set_xticks(np.arange(0, len_ys, 1))
                plt.yticks(np.arange(0, 100, 1))
                ax.plot(xs_normal, ys_normal, color="r", label="normal")
                ax.plot(xs_custom, ys_custom, color="b", label="custom")
                ax.legend()
                plt.grid(axis="y", linestyle="-")
                plt.savefig("stats_"+s+".png")
                plt.close()
        
def main(filesToProcess, posiciones):
    print(posiciones)
    data_positions = []
    data_interval = []
    if(("-" in posiciones)):
        data_interval = [int(x) for x in posiciones.split("-")]
        pos_init = data_interval[0]
        pos_end = data_interval[1]
        for pos in range (pos_init,pos_end+1):
            data_positions.append(pos)
    else:
        data_positions = [int(x) for x in posiciones.split(",")]
    record = Record()
    cnt = 1
    for filepath in filesToProcess:
        try:
            with open(filepath, 'r') as fp:
                
                line = fp.readline()
                
                while line:
                    line = fp.readline()
                    stage = record.current_stage
                    record.check_stage(line)
                    # solo para comprobar que esta encontrando las fases
                    if (stage != record.current_stage):
                        print("change stage to ",record.current_stage)
                    
                    record.record_stat(line)
                    cnt +=1

                print("log lines", cnt)
                record.print_stats()
                #record.save_bar_plots()
                #record.save_data_plots()
        finally:
            fp.close()
    record.save_bar_plots(data_positions)

    
    
if __name__ == "__main__":
    files = ["salida_1.txt","salida_2.txt", "salida_3.txt", "salida_4.txt"]
    main(files, sys.argv[1])


