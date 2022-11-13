ifneq ($(OPT),opt-gccold)
ifneq ($(OPT),dbg)
OFLAGS          +=      -mtune=native
endif
endif
CFLAGS		+=	-m32
CXXFLAGS	+=	-m32
LDFLAGS		+=      -m32

ifeq ($(SHARED),true)
LIBBUILDFLAGS	+=     	-m32
endif
