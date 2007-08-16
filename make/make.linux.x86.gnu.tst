#--- $Id: make.linux.x86.gnu.tst,v 1.1 2007/08/16 13:47:38 bzfpfend Exp $
FLAGS		+=	-DNDEBUG # -DROUNDING_FE
OFLAGS		+=	-g -O0 -march=pentiumpro -fomit-frame-pointer # -malign-double -mcpu=pentium4 -g
CFLAGS		+=	$(GCCWARN) -Wno-strict-aliasing -Wno-missing-declarations -Wno-missing-prototypes
CXXFLAGS	+=	$(GXXWARN) -Wno-strict-aliasing # -fno-exceptions (CLP uses exceptions)
ARFLAGS		=	crs
LDFLAGS		+=
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses
