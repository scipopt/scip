ARFLAGS		=	crs
LDFLAGS		+=
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp
MPFR_LDFLAGS	=	-lmpfr
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses

ifeq ($(LPS),cpx)
LPSLDFLAGS	+=	  -Wl,-no_compact_unwind -framework IOKit -framework Carbon
endif

ifeq ($(SHARED),true)
LIBBUILDFLAGS   =	-dynamiclib -undefined suppress -flat_namespace
FLAGS		+=	-fPIC
endif
LINK_shared		=	-shared

CXXFLAGS	=	-ffp-contract=off -fno-stack-check
CFLAGS		=	-ffp-contract=off -std=c99 -fno-stack-check -D_XOPEN_SOURCE=600

ifeq ($(LTO),true)
  # for GCC < 10, use just -flto, otherwise use -flto=auto, which should give faster link times
  # -fno-fat-lto-objects (since GCC 5) should improve compilation time a bit
  GCCVERSION := $(shell $(CC) -dumpversion | cut -f1 -d.)
  LTOFLAG := $(word $(shell expr \( $(GCCVERSION) \>= 10 \) + 1), -flto -flto=auto)

  CFLAGS	+=	$(LTOFLAG) -fno-fat-lto-objects
  CXXFLAGS	+=	$(LTOFLAG) -fno-fat-lto-objects
  LDFLAGS	+=	$(LTOFLAG) -Wno-stringop-overflow -Wno-alloc-size-larger-than -Wno-odr
ifeq ($(SHARED),true)
  LIBBUILDFLAGS +=	$(LTOFLAG) -Wno-stringop-overflow -Wno-alloc-size-larger-than -Wno-odr
endif
endif

FORTRANLIBS	=	-lgfortran
FORTRAN_NAMING_CONVENTION = LCASE_DECOR
