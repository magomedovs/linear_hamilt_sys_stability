import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

color_dict = {'stable' : 'green', 'unstable' : 'red', 
	'double_minus_one' : 'blue', 'double_plus_one' : 'black'}

filename = sys.argv[1] 
if (os.path.exists(filename)):
	df = pd.read_table(filename, delimiter=r"\s+")
	fig, ax = plt.subplots(nrows=1, ncols=1)
	
	colors = [color_dict[stab_stat] for stab_stat in df['stability_status']]
	ax.scatter(df['beta'], df['h'], c=colors, linewidths=0.1)
	
	ax.grid(True)
	plt.show()

else:
	print("File not found")

