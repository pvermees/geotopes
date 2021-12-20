fnames <- list.files(pattern = "\\.R$")
for (fname in fnames){
    if (fname!='test.R') source(fname)
}
