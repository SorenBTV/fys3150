import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def histogram(filename):

    data = pd.read_csv(filename, sep=" ", header=1)

    return print(data())


histogram("problem_5_T1.0_cycles10000_not_ordered.txt")
