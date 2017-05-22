ARFLAGS		=	crs
LDFLAGS		+=
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses

ifeq ($(LPS),cpx)
LPSLDFLAGS	+=	  -Wl,-no_compact_unwind -framework IOKit -framework Carbon
endif

ifeq ($(SHARED),true)
LIBBUILDFLAGS   =	-dynamiclib -undefined suppress -flat_namespace
endif

CXXFLAGS	=	-std=c++0x
CFLAGS		=	-std=c99 -D_XOPEN_SOURCE=600

FORTRANLIBS	=	-lgfortran
FORTRAN_NAMING_CONVENTION = LCASE_DECOR
