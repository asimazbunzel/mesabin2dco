# This is the makefile for the package template

#################################################################
#
# PREAMBLE

include $(MESA_DIR)/utils/makefile_header

#################################################################
#
# SOURCE FILES

SRCS = \
	cc_def.f90 \
	cc_data.f90 \
	cc_mechanism.f90 \
	cc_ctrls_io.f90 \
	cc_support.f90 \
	cc_lib.f90

#################################################################
#
# LIBRARIES

ifeq ($(USE_SHARED), YES)
	LIBS_OTHER = $(LIBS_MESA_STAR)
	DEPS_OTHER = $(patsubst %,$(MESA_LIB_DIR)/lib%.$(LIB_SUFFIX),$(LIBS_OTHER))
	LOAD_OTHER = -L$(MESA_LIB_DIR) $(LOAD_MESA_STAR)
endif

#################################################################
#
# TARGETS

LIB = libcc.a

OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))

$(LIB) : $(OBJS) $(DEPS_OTHER)
ifneq ($(QUIET),)
	@echo LIB_TOOL $(LIB)
	@$(LIB_TOOL) $(LIB) $(OBJS) $(LOAD_OTHER)
else
	$(LIB_TOOL) $(LIB) $(OBJS) $(LOAD_OTHER)
endif

clean:
	-@rm -f *.o *.mod *.modmic *genmod.f90 *.so *.a .depend *.smod

MKDIR_P = mkdir -p

install:
	@(MKDIR_P) ../../../include ../../../lib
	@$(CP_IF_NEWER) cc_lib.mod ../../../include
	@$(CP_IF_NEWER) cc_def.mod ../../../include
	@$(CP_IF_NEWER) $(LIB) ../../../lib

nodeps : $(.DEFAULT_GOAL)

#################################################################
#
# COMPILATION RULES

COMPILE = $(COMPILE_TO_DEPLOY) $(FCfree)
#COMPILE = $(COMPILE_TO_TEST) $(FCfree)

COMPILE_CMD = $(COMPILE)

%.o : %.mod

%.o : %.f
ifneq ($(QUIET),)
	@echo COMPILE_CMD $<
	@$(COMPILE_CMD) $<
else
	$(COMPILE_CMD) $<
endif

%.o : %.f90
ifneq ($(QUIET),)
	@echo COMPILE_CMD $<
	@$(COMPILE_CMD) $<
else
	$(COMPILE_CMD) $<
endif

%.mod : %.o
	@true

#################################################################
#
# DEPENDENCIES

SRC_PATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR):$(MOD_OTHER_DIR)

vpath %.f90 $(SRC_PATH)
vpath %.inc $(SRC_PATH)

vpath %.mod $(MESA_DIR)/include

NODEPS = $(or $(filter nodeps,$(MAKECMDGOALS)),$(filter clean,$(MAKECMDGOALS)))

ifeq ($(NODEPS),)

  .depend :
  ifneq ($(QUIET),)
	@echo MAKEDEPF90
	@$(MAKEDEPF90) -I$(SRC_PATH) $(SRCS) > .depend
  else
	$(MAKEDEPF90) -I$(SRC_PATH) $(SRCS) > .depend
  endif

  -include .depend

endif
