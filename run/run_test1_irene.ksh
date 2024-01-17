#!/bin/bash

ln -sf ../src/qg.e .
ln -sf ../test/params_01.in params.in

chmod +x job_test1_irene.ksh
ccc_msub job_test1_irene.ksh
