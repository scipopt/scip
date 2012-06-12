ARFLAGS		=	crs
LDFLAGS		+=
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp -lgmpxx
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses

ifeq ($(LPS),cpx)
LPSLDFLAGS	+=	 -framework IOKit -framework Carbon
endif

ifeq ($(SHARED),true)
LIBBUILDFLAGS   =	-dynamiclib -undefined suppress -flat_namespace
endif