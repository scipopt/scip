ARFLAGS		=	crs
ifeq ($(IPOPT),true)
LDFLAGS                +=      -Wl,--no-as-needed
endif
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp
MPFR_LDFLAGS	=	-lmpfr
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses
CFLAGS		=	-ffp-contract=off -std=c99 -D_XOPEN_SOURCE=600
CXXFLAGS	=	-ffp-contract=off
FORTRANLIBS	=	-lgfortran
FORTRAN_NAMING_CONVENTION = LCASE_DECOR
ifeq ($(SHARED),true)
FLAGS		+=	-fPIC
endif
LINK_shared		=	-shared

# ThreadSanitizer (https://github.com/google/sanitizers/wiki/ThreadSanitizerCppManual)
ifeq ($(SANITIZE),thread)
  SANITIZERFLAGS = -g -fsanitize=thread
endif

# AddressSanitizer (https://github.com/google/sanitizers/wiki/AddressSanitizer)
ifeq ($(SANITIZE),address)
  SANITIZERFLAGS = -g -fsanitize=address
endif

ifeq ($(SANITIZE),memory)
  $(warning Memory Sanitizer not available with GCC)
endif

# UndefinedBehaviorSanitizer if SANITIZE is true, thread, address, or memory
ifneq ($(filter $(SANITIZE),true thread address memory),)
  SANITIZERFLAGS += -g -fsanitize=undefined -fsanitize=float-cast-overflow -fsanitize=float-divide-by-zero
endif

CFLAGS += $(SANITIZERFLAGS)
CXXFLAGS += $(SANITIZERFLAGS)
LDFLAGS += $(SANITIZERFLAGS)
ifeq ($(SHARED),true)
  LIBBUILDFLAGS += $(SANITIZERFLAGS)
endif

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
