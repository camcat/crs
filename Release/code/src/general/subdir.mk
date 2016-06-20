################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/general/CRS_LogLikelihood.c \
../code/src/general/calculateDCFSperturbed.c \
../code/src/general/eqkfm_copy.c \
../code/src/general/find_timesteps.c \
../code/src/general/forecast_stepG.c \
../code/src/general/lin_interp_eqkfm.c \
../code/src/general/mem_mgmt.c \
../code/src/general/setup.c \
../code/src/general/setup_time.c \
../code/src/general/struct_conversions.c 

OBJS += \
./code/src/general/CRS_LogLikelihood.o \
./code/src/general/calculateDCFSperturbed.o \
./code/src/general/eqkfm_copy.o \
./code/src/general/find_timesteps.o \
./code/src/general/forecast_stepG.o \
./code/src/general/lin_interp_eqkfm.o \
./code/src/general/mem_mgmt.o \
./code/src/general/setup.o \
./code/src/general/setup_time.o \
./code/src/general/struct_conversions.o 

C_DEPS += \
./code/src/general/CRS_LogLikelihood.d \
./code/src/general/calculateDCFSperturbed.d \
./code/src/general/eqkfm_copy.d \
./code/src/general/find_timesteps.d \
./code/src/general/forecast_stepG.d \
./code/src/general/lin_interp_eqkfm.d \
./code/src/general/mem_mgmt.d \
./code/src/general/setup.d \
./code/src/general/setup_time.d \
./code/src/general/struct_conversions.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/general/%.o: ../code/src/general/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc   -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


