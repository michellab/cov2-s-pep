#!/usr/bin/python3

import csv
import glob
import numpy as np
from datetime import datetime
fmt="%d %H:%M:%S"
sim_dirs = glob.glob("HR*")

#iterate through folders:
for sim_folder in sim_dirs:
    
    steps = []
    write_moments = []
    with open(sim_folder+"/md.log") as md_logfile:
        reader = csv.reader(md_logfile)
        for line in reader:
            try:
                if line[0].startswith("Writing"):
                    step = int(line[1].rsplit(" ")[2])
                    steps.append(step)
                    write_moment = line[1].rsplit(" ")[6]+" "+line[1].rsplit(" ")[7]
                    write_moments.append(write_moment)
            except:
                pass
    # now compute elapsed time per step compared to step 0:
    step0_date = write_moments[0]
    elapsed_times = []
    for write_moment in write_moments:
        tstamp0 = datetime.strptime(step0_date, fmt)
        tstamp1 = datetime.strptime(write_moment, fmt)
        td = tstamp1 - tstamp0
        hours_passed =  td.total_seconds() / 60 / 60
        elapsed_times.append(hours_passed)
    
    # based on step vs time elapsed, compute predicted total simulation time per step and compute mean:
    total_hours_predictions = []
    for step, timestamp in zip(steps, elapsed_times):
        predicted_total_hours = 50000000/step * timestamp
        total_hours_predictions.append(predicted_total_hours)
    mean_predicted_runtime = np.mean(total_hours_predictions[10:])


    #subtract from time already elapsed:
    time_left = mean_predicted_runtime - elapsed_times[-1]
    print(sim_folder, "finishes in", round(time_left, 2), "hours." )
