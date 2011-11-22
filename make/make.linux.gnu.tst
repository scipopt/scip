FLAGS		+=	-DNDEBUG -DROUNDING_FE
OFLAGS		+=	-g -O0 -fomit-frame-pointer # -malign-double -mcpu=pentium4 -g
CFLAGS		+=	$(GCCWARN) -Wno-strict-aliasing -Wno-missing-declarations -Wno-missing-prototypes
CXXFLAGS	+=	$(GXXWARN) -Wno-strict-aliasing # -fno-exceptions (CLP uses exceptions)
