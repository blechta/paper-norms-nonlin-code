#!/bin/bash

python main.py NicaiseVenel -0.33333333333333333 08 | tee NicaiseVenel_-0.333_08.log
python main.py NicaiseVenel -0.33333333333333333 16 | tee NicaiseVenel_-0.333_16.log
python main.py NicaiseVenel -0.33333333333333333 32 | tee NicaiseVenel_-0.333_32.log

python main.py NicaiseVenel -0.01 08 | tee NicaiseVenel_-0.01_08.log
python main.py NicaiseVenel -0.01 16 | tee NicaiseVenel_-0.01_16.log
python main.py NicaiseVenel -0.01 32 | tee NicaiseVenel_-0.01_32.log

python main.py NicaiseVenel -0.99 08 | tee NicaiseVenel_-0.99_08.log
python main.py NicaiseVenel -0.99 16 | tee NicaiseVenel_-0.99_16.log
python main.py NicaiseVenel -0.99 32 | tee NicaiseVenel_-0.99_32.log


grep ^RESULT *.log          > out.txt
grep -m 1 Estimators *.log >> out.txt
grep "nabla r" *.log       >> out.txt
