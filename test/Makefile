testall = test_2D test_3D

# skip 3D tests in CI
ifeq ($(TRAVIS),true)
	testall = test_2D
endif

all :: $(testall)

test_2D ::
	@echo "Testing all in 2D ..."
	@cd 2D; make

test_3D ::
	@echo "Testing all in 3D ..."
	@cd 3D; make
