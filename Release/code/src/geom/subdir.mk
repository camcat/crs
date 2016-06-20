################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/geom/check_same_geometry.c \
../code/src/geom/convert_geometry.c \
../code/src/geom/coord_trafos.c \
../code/src/geom/find_gridpoints.c \
../code/src/geom/top_of_slipmodel.c 

OBJS += \
./code/src/geom/check_same_geometry.o \
./code/src/geom/convert_geometry.o \
./code/src/geom/coord_trafos.o \
./code/src/geom/find_gridpoints.o \
./code/src/geom/top_of_slipmodel.o 

C_DEPS += \
./code/src/geom/check_same_geometry.d \
./code/src/geom/convert_geometry.d \
./code/src/geom/coord_trafos.d \
./code/src/geom/find_gridpoints.d \
./code/src/geom/top_of_slipmodel.d 


# Each subdirectory must supply rules for building sources it contributes
code/src/geom/%.o: ../code/src/geom/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc   -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


