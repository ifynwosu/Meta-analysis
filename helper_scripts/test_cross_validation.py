
# import relevant libraries

import pandas as pd
from sklearn import datasets, linear_model
from sklearn.model_selection import train_test_split
# from matplotlib import pyplot as plt
import numpy as np



indicies = []
train_test_set = []
last_element = 0
for j in range(len(datasets)):
    train_test_set.append(datasets[j])
    indicies.append(np.arange(last_element,last_element+len(datasets[j])))
    last_element += len(datasets[j])

cv_list = []
for i in range(15):
    cv_train = np.hstack([indicies[x] for x in range(15) if x != i])
    cv_list.append((cv_train,indicies[i]))


# wget https://www.openssl.org/source/openssl-1.1.1u.tar.gz
# tar xzf openssl-1.1.1u.tar.gz
# cd openssl-1.1.1u
# ./config --prefix=/home/inwosu/packages/openssl --openssldir=/home/inwosu/packages/openssl shared -Wl,--enable-new-dtags,-rpath,'$(LIBRPATH)'


# wget https://www.python.org/ftp/python/3.11.4/Python-3.11.4.tgz
# tar xf Python-3.9.2.tgz
# cd Python-3.11.4
# ./configure --prefix=/home/inwosu/packages/python311 --enable-optimizations --with-openssl=/home/inwosu/packages/openssl

# CFLAGS="-I/home/inwosu/packages/openssl/include/" LDFLAGS="${LDFLAGS} -Wl,-rpath=/home/inwosu/packages/openssl/lib" ./configure --prefix=/home/inwosu/packages/python311 --enable-optimizations --with-openssl=/home/inwosu/packages/openssl
# make
# make install 


# ldd build/lib.linux-x86_64-3.9/_ssl.cpython-39-x86_64-linux-gnu.so  | grep weird


# WARNING: The scripts pip3 and pip3.11 are installed in '/home/inwosu/packages/python311/bin' which is not on PATH.
# Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
# Successfully installed pip-23.1.2 setuptools-65.5.0

# /usr/bin/python3
