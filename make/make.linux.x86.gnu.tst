#--- $Id: make.linux.x86.gnu.tst,v 1.5 2008/09/30 18:00:03 bzfberth Exp $
FLAGS		+=	-DNDEBUG # -DROUNDING_FE
OFLAGS		+=	-g -O0 -march=pentiumpro -fomit-frame-pointer # -malign-double -mcpu=pentium4 -g
CFLAGS		+=	$(GCCWARN) -Wno-strict-aliasing -Wno-missing-declarations -Wno-missing-prototypes
CXXFLAGS	+=	$(GXXWARN) -Wno-strict-aliasing # -fno-exceptions (CLP uses exceptions)
ARFLAGS		=	crs
LDFLAGS		+=	-Wl,-rpath,$(CURDIR)/$(LIBDIR)
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses
