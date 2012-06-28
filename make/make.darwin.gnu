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