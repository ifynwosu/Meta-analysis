import os, sys

# Data Processing
import pandas as pd
import numpy as np

from numpy import mean
from numpy import std

# scikit learn modules
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict