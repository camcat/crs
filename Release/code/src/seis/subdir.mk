################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/seis/GR.c \
../code/src/seis/WellsCoppersmith.c \
../code/src/seis/background_rate.c \
../code/src/seis/cmbopt.c \
../code/src/seis/decluster.c \
../code/src/seis/smoothed_rate_Helmstetter.c \
../code/src/seis/soumod1.c 

OBJS += \
./code/src/seis/GR.o \
./code/src/seis/WellsCoppersmith.o \
./code/src/seis/background_rate.o \
./code/src/seis/cmbopt.o \
./code/src/seis/decluster.o \
./code/src/seis/smoothed_rate_Helmstetter.o \
./code/src/seis/soumod1.o 

C_DEPS += \
./code/src/seis/GR.d \
./code/src/seis/WellsCoppersmith.d \
./code/src/seis/background_rate.d \
./code/src/seis/cmbopt.d \
./code/src/seis/decluster.d \
./code/src/seis/smoothed_rate_Helmstetter.d \
./code/src/seis/soumod1.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/seis/%.o: ../code/src/seis/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc   -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


