import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_table("problem_5_T1.0_cycles10000_not_ordered.txt", header = None)
data2 = pd.read_table("problem_5_T1.0_cycles10000_ordered.txt", header = None)
data3 = pd.read_table("problem_5_T2.4_cycles10000_not_ordered.txt", header = None)
data4 = pd.read_table("problem_5_T2.4_cycles10000_ordered.txt", header = None)

e = []
E = []
e2 = []
E2 = []



for i in range(0, len(data[0])):
    a = data[0][i].split()
    e.append(float(a[1]))
    E.append(float(a[3]/))

for i in range(0, len(data3[0])):
    b = data3[0][i].split()
    e2.append(float(b[1]))

plt.hist(e, bins=10, label="T=1.0")
plt.hist(e2, bins=10, label ="T=2.4")
plt.legend()

plt.show()
