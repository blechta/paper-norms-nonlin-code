#!/bin/bash

python main.py NicaiseVenel -0.33333333333333333 008 | tee NicaiseVenel_-0.333_008.log
python main.py NicaiseVenel -0.33333333333333333 016 | tee NicaiseVenel_-0.333_016.log
python main.py NicaiseVenel -0.33333333333333333 032 | tee NicaiseVenel_-0.333_032.log
python main.py NicaiseVenel -0.33333333333333333 064 | tee NicaiseVenel_-0.333_064.log
python main.py NicaiseVenel -0.33333333333333333 128 | tee NicaiseVenel_-0.333_128.log

python main.py NicaiseVenel -0.01 008 | tee NicaiseVenel_-0.01_008.log
python main.py NicaiseVenel -0.01 016 | tee NicaiseVenel_-0.01_016.log
python main.py NicaiseVenel -0.01 032 | tee NicaiseVenel_-0.01_032.log
python main.py NicaiseVenel -0.01 064 | tee NicaiseVenel_-0.01_064.log
python main.py NicaiseVenel -0.01 128 | tee NicaiseVenel_-0.01_128.log

python main.py NicaiseVenel -0.99 008 | tee NicaiseVenel_-0.99_008.log
python main.py NicaiseVenel -0.99 016 | tee NicaiseVenel_-0.99_016.log
python main.py NicaiseVenel -0.99 032 | tee NicaiseVenel_-0.99_032.log
python main.py NicaiseVenel -0.99 064 | tee NicaiseVenel_-0.99_064.log
python main.py NicaiseVenel -0.99 128 | tee NicaiseVenel_-0.99_128.log

python main.py BonnetBenDhia -5 008 | tee BonnetBenDhia_-5_008.log
python main.py BonnetBenDhia -5 016 | tee BonnetBenDhia_-5_016.log
python main.py BonnetBenDhia -5 032 | tee BonnetBenDhia_-5_032.log
python main.py BonnetBenDhia -5 064 | tee BonnetBenDhia_-5_064.log
python main.py BonnetBenDhia -5 128 | tee BonnetBenDhia_-5_128.log

grep ^RESULT *.log          > out.txt
grep -m 1 Estimators *.log >> out.txt
grep "nabla r" *.log       >> out.txt
grep "nabla(u-u_h)" *.log  >> out.txt
