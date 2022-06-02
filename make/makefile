ifeq ($(MESA_DIR),)
ifeq ($($MESA_DIR_INTENTIONALLY_EMPTY),)
  $(error MESA_DIR enviroment variable is not set)
endif
endif

include $(MESA_DIR)/utils/makefile_header

# we do not use MESA default include and lib folders
# instead, we have local ones with our modules cc and ce
OWN_INCLUDE_DIR = ../include
OWN_LIB_DIR = ../lib

OWN_LIBS = cc ce $(LIBS_MESA_BINARY)
LOAD_OWN = -L$(MESA_LIB_DIR) -L$(OWN_LIB_DIR) $(addprefix -l,$(OWN_LIBS))
LOAD_ALL = $(LOAD_OWN) $(LOAD_MESA_BINARY_EXT)

#################################################################

OWN_WORK_COMPILE = \
	$(FC) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCdebug) $(FCfree) \
	-I$(MESA_INCLUDE_DIR) -I$(OWN_INCLUDE_DIR) $(INCLUDE_HDF5) -c

ifndef STAR
STAR = binary
endif

OBJS = $(BINARY_WORK_OBJS) binary_run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job
BINARY_JOB_DIR = $(MESA_DIR)/binary/job

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_ALL)
	
#################################################################

# WORK_COMPILE defined in makefile_header

# STAR_WORK_OBJS

run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(OWN_WORK_COMPILE) $<

run_star.o: $(STAR_JOB_DIR)/run_star.f90
	$(OWN_WORK_COMPILE) $<


# BINARY_WORK_OBJS

binary_run.o: $(WORK_SRC_DIR)/binary_run.f90
	$(OWN_WORK_COMPILE) $<

run_binary_extras.o: $(WORK_SRC_DIR)/run_binary_extras.f90
	$(OWN_WORK_COMPILE) $<

run_binary.o: $(BINARY_JOB_DIR)/run_binary.f90
	$(OWN_WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod *.smod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
