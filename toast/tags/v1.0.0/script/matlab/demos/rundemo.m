function rundemo(demo)
tgt = which(demo);
[pathstr name] = fileparts(tgt);
eval(['cd ' ['''' pathstr '''']]);
eval(name)