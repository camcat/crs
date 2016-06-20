################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../code/src/inp_out/print_output.c \
../code/src/inp_out/read_crust.c \
../code/src/inp_out/read_csep_template.c \
../code/src/inp_out/read_eqkfm.c \
../code/src/inp_out/read_eqkfm_fsp.c \
../code/src/inp_out/read_focmec.c \
../code/src/inp_out/read_inputfile.c \
../code/src/inp_out/read_matrix.c \
../code/src/inp_out/read_param.c \
../code/src/inp_out/read_zmap.c \
../code/src/inp_out/write_csep_forecast.c

OBJS += \
./code/src/inp_out/print_output.o \
./code/src/inp_out/read_crust.o \
./code/src/inp_out/read_csep_template.o \
./code/src/inp_out/read_eqkfm.o \
./code/src/inp_out/read_eqkfm_fsp.o \
./code/src/inp_out/read_focmec.o \
./code/src/inp_out/read_inputfile.o \
./code/src/inp_out/read_matrix.o \
./code/src/inp_out/read_param.o \
./code/src/inp_out/read_zmap.o \
./code/src/inp_out/write_csep_forecast.o

C_DEPS += \
./code/src/inp_out/print_output.d \
./code/src/inp_out/read_crust.d \
./code/src/inp_out/read_csep_template.d \
./code/src/inp_out/read_eqkfm.d \
./code/src/inp_out/read_eqkfm_fsp.d \
./code/src/inp_out/read_focmec.d \
./code/src/inp_out/read_inputfile.d \
./code/src/inp_out/read_matrix.d \
./code/src/inp_out/read_param.d \
./code/src/inp_out/read_zmap.d \
./code/src/inp_out/write_csep_forecast.d


# Each subdirectory must supply rules for building sources it contributes
code/src/inp_out/%.o: ../code/src/inp_out/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc   -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


