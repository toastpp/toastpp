MATLAB = matlab -nodesktop -nojvm

all :: fwd1 test1 test2

fwd1 ::
	@echo -----------------------------------------------------------
	@echo matlab forward solver test 1 of 1 ...
	@echo Comparison Explicit/direct moments
	@echo -----------------------------------------------------------
	$(MATLAB) -r "istest=1;fwd_tpsf;if pass,exit,end"
	@echo "#### Test 1 passed! ####"


test1 ::
	@echo -----------------------------------------------------------
	@echo matlab toast test 1 of 2 ...
	@echo Gauss-Newton solver
	@echo -----------------------------------------------------------
	$(MATLAB) -r "recon1();exit"
	cmp recon.dat recon1.dat
	@rm recon.dat
	@echo "#### Test 1 passed! ####"

test2 ::
	@echo -----------------------------------------------------------
	@echo matlab toast test 2 of 2 ...
	@echo conjugate gradient solver
	@echo -----------------------------------------------------------
	$(MATLAB) -r "recon2();exit"
	cmp recon.dat recon2.dat
	@rm recon.dat
	@echo "#### Test 2 passed! ####"
