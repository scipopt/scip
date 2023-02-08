ARFLAGS		=	crs
ifeq ($(IPOPT),true)
LDFLAGS                +=      -Wl,--no-as-needed
endif
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses
CFLAGS		=	-std=c99 -D_XOPEN_SOURCE=600
FORTRANLIBS	=	-lgfortran
FORTRAN_NAMING_CONVENTION = LCASE_DECOR
ifeq ($(SHARED),true)
FLAGS		+=	-fPIC
endif
LINK_shared		=	-shared
