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

python main.py BonnetBenDhia -3.1 008 | tee BonnetBenDhia_-3.1_008.log
python main.py BonnetBenDhia -3.1 016 | tee BonnetBenDhia_-3.1_016.log
python main.py BonnetBenDhia -3.1 032 | tee BonnetBenDhia_-3.1_032.log
python main.py BonnetBenDhia -3.1 064 | tee BonnetBenDhia_-3.1_064.log
python main.py BonnetBenDhia -3.1 128 | tee BonnetBenDhia_-3.1_128.log

python main.py BonnetBenDhia_adaptive -5   008 | tee BonnetBenDhia_adaptive_-5_008.log
python main.py BonnetBenDhia_adaptive -3.1 008 | tee BonnetBenDhia_adaptive_-3.1_008.log


grep -h ^RESULT *.log > out.txt
