function smat = toastDotSysmat (hmesh,mua,mus,ref,freq)
smat = toast(uint32(23),hmesh.handle,double(mua),double(mus),double(ref),freq);
end
