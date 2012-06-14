ARFLAGS		=	crs
ifeq ($(IPOPT),true)
LDFLAGS		+=      -Wl,--no-as-needed
endif
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp -lgmpxx
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses
