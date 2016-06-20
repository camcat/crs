################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/okada/dc3d.c \
../code/src/okada/okadaDCFS.c \
../code/src/okada/prestress.c \
../code/src/okada/pscokada.c 

OBJS += \
./code/src/okada/dc3d.o \
./code/src/okada/okadaDCFS.o \
./code/src/okada/prestress.o \
./code/src/okada/pscokada.o 

C_DEPS += \
./code/src/okada/dc3d.d \
./code/src/okada/okadaDCFS.d \
./code/src/okada/prestress.d \
./code/src/okada/pscokada.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/okada/%.o: ../code/src/okada/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc   -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


