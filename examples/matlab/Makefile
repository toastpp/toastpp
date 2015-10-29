matlab_doc :: matlab_doc_fwd1 matlab_doc_recon_cg

matlab_doc_fwd1 ::
	rm -rf html
	matlab -nodesktop -r "publish('fwd1.m');exit"
	mv html $(TOASTDIR)/script/matlab/html/fwd1

matlab_doc_recon_cg ::
	rm -rf html
	matlab -nodesktop -r "publish('recon_cg.m');exit"
	mv html $(TOASTDIR)/script/matlab/html/recon_cg
