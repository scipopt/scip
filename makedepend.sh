# disable ZIMPL for dependency generation
make LPS=cpx OPT=opt ZIMPL=false depend
make LPS=cpx OPT=dbg ZIMPL=false depend
make LPS=spx OPT=opt ZIMPL=false lpidepend
make LPS=spx OPT=dbg ZIMPL=false lpidepend
make LPS=spx121 OPT=opt ZIMPL=false lpidepend
make LPS=spx121 OPT=dbg ZIMPL=false lpidepend
make LPS=clp OPT=opt ZIMPL=false lpidepend
make LPS=clp OPT=dbg ZIMPL=false lpidepend
#make LPS=xprs OPT=opt lpidepend
#make LPS=xprs OPT=dbg lpidepend
