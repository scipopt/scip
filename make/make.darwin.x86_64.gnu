ifneq ($(OPT),opt-gccold)
ifneq ($(OPT),dbg)
OFLAGS          +=      -mtune=native  # -malign-double -mcpu=pentium4
endif
endif
CFLAGS		+=	-m64  -fno-stack-check
CXXFLAGS	+=	-m64  -fno-stack-check
LDFLAGS		+=      -m64

ifeq ($(SHARED),true)
LIBBUILDFLAGS   +=	-m64
endif
