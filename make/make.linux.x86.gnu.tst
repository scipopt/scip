FLAGS		+=	-DNDEBUG -DROUNDING_FE
OFLAGS		+=	-g -O0 -fomit-frame-pointer # -malign-double -mcpu=pentium4 -g
CFLAGS		+=	-m32 $(GCCWARN) -Wno-strict-aliasing -Wno-missing-declarations -Wno-missing-prototypes
CXXFLAGS	+=	-m32 $(GXXWARN) -Wno-strict-aliasing # -fno-exceptions (CLP uses exceptions)
ARFLAGS		=	crs
LDFLAGS		+=	-m32
ZLIB_FLAGS	=
ZLIB_LDFLAGS 	=	-lz
GMP_FLAGS	=
GMP_LDFLAGS 	=	-lgmp
READLINE_FLAGS	=
READLINE_LDFLAGS=	-lreadline -lncurses
