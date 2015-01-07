SHELL = /bin/bash
MATLAB = matlab -nodesktop -nojvm

all :: fwd1 test1 test2

fwd1 ::
	@echo -----------------------------------------------------------
	@echo matlab forward solver test 1 of 1 ...
	@echo Comparison Explicit/direct moments
	$(MATLAB) -r "run([getenv('TOASTDIR') '/mtoast_install.m']);istest=1;fwd_tpsf;if pass,exit,end"
	@echo "#### Test 1 passed! ####"
	@echo -----------------------------------------------------------


test1 ::
	@echo -----------------------------------------------------------
	@echo matlab toast test 1 of 2 ...
	@echo Gauss-Newton solver
	$(MATLAB) -r "recon1" &> /dev/null
	@echo "#### Test 1 passed! ####"
	@echo -----------------------------------------------------------

test2 ::
	@echo -----------------------------------------------------------
	@echo matlab toast test 2 of 2 ...
	@echo conjugate gradient solver
	$(MATLAB) -r "recon2();exit"
	cmp recon.dat recon2.dat
	@rm recon.dat
	@echo "#### Test 2 passed! ####"
	@echo -----------------------------------------------------------
