################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include code/src/nr/subdir.mk
-include code/src/util/subdir.mk
-include code/src/seis/subdir.mk
-include code/src/okada/subdir.mk
-include code/src/inp_out/subdir.mk
-include code/src/geom/subdir.mk
-include code/src/general/subdir.mk
-include code/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: run_crs

# Tool invocations
run_crs: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -o "run_crs" $(OBJS) $(USER_OBJS) $(LIBS) -lm -lgomp -lgsl -lgslcblas
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) run_crs
	-$(RM) code/src/nr/*o code/src/nr/*d
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
