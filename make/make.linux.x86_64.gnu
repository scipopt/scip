ifneq ($(OPT),opt-gccold)
ifneq ($(OPT),dbg)
OFLAGS          +=      -mtune=native
endif
endif
CFLAGS		+=	-m64
CXXFLAGS	+=	-m64
LDFLAGS		+=      -m64

ifeq ($(SHARED),true)
LIBBUILDFLAGS	+=     	-m64
endif
