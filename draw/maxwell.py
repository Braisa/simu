import numpy as np
import matplotlib.pyplot as plt
import pandas as pd # type: ignore
import os

path = os.path.dirname(__file__) + "/"

rva = pd.read_csv(path + "big_rva.txt", sep=" ")

rva["v2"] = rva["vx"]**2 + rva["vy"]**2 + rva["vz"]**2
