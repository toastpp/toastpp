include $(TOASTDIR)/Makefile.incl

DGDIR = /cs/research/medim/projects2/projects/sprerapa/trunk/src/libdgfe
DGINC = -I$(DGDIR) -I$(DGDIR)/../libfe 
DGLIB = -L$(TLIB) -ldgfe
TOASTLIB = -L$(TLIB) -lstd -lstoast -lsuperlu -lfe -lblzpack -larpack -llapack $(BLASLIB) $(SPBLASLIB) $(FLIBS)
BIN = .

.cc.o :
	$(CXX) $(STOASTFLAGS) $(DGINC) -c $<

test : 	test.cc nonconformingMesh.cc dgfwdsolver.cc dgSource.cc
	g++ -I/cs/research/medim/projects2/projects/sprerapa/trunk -I/cs/research/medim/projects2/projects/sprerapa/trunk/include -I/cs/research/medim/projects2/projects/sprerapa/trunk/src/libmath -I/cs/research/medim/projects2/projects/sprerapa/trunk/src/libfe -I/cs/research/medim/projects2/projects/sprerapa/trunk/src/libstoast test.cc nonconformingMesh.cc dgfwdsolver.cc dgSource.cc -L/cs/research/medim/projects2/projects/sprerapa/trunk/gcc4.1.2_shared/lib -lfe -lmath -lsuperlu -lblzpack -larpack -llapack -lblas -lfsbtk -lstoast -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../../../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/../../.. -L/lib/../lib64 -L/usr/lib/../lib64 -lg2c -o test
	rm -rf *.o

mytest: test.cc nonconformingMesh.cc
	g++ $(DGINC) $(DGLIB) $^ -o $(BIN)/$@

clean:: 
	$(RM) *.o

