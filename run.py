import sys
import os
import pandas as pd
import json

import matplotlib.pyplot as plt
from calculator import calculator
paths = os.getcwd() + "/" 
sys.path.insert(0, paths + "/src")
from EDA import eda

def main(targets):
  



    if ("test" in targets) or (targets == "test):
        print(Given )
if __name__ == '__main__':
    
    targets = sys.argv[1:]
    
    main(targets)
