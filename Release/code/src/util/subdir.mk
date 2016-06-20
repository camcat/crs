################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/util/util1.c \
../code/src/util/error.c \
../code/src/util/files.c \
../code/src/util/fit_splines.c \
../code/src/util/interp_quad.c \
../code/src/util/merge.c \
../code/src/util/moreutil.c \
../code/src/util/mscorr.c \
../code/src/util/roots3.c \
../code/src/util/splines_eqkfm.c 

OBJS += \
./code/src/util/util1.o \
./code/src/util/error.o \
./code/src/util/files.o \
./code/src/util/fit_splines.o \
./code/src/util/interp_quad.o \
./code/src/util/merge.o \
./code/src/util/moreutil.o \
./code/src/util/mscorr.o \
./code/src/util/roots3.o \
./code/src/util/splines_eqkfm.o 

C_DEPS += \
./code/src/util/util1.d \
./code/src/util/error.d \
./code/src/util/files.d \
./code/src/util/fit_splines.d \
./code/src/util/interp_quad.d \
./code/src/util/merge.d \
./code/src/util/moreutil.d \
./code/src/util/mscorr.d \
./code/src/util/roots3.d \
./code/src/util/splines_eqkfm.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/util/%.o: ../code/src/util/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -D_no_numericalrecipes -I/usr/local/include -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


