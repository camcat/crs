################################################################################
# Automatically-generated file. Do not edit!
################################################################################

ifdef numrec_path

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
$(numrec_path)/recipes/gasdev.c \
$(numrec_path)/recipes/gaussj.c \
$(numrec_path)/recipes/spline.c \
$(numrec_path)/recipes/splint.c \
$(numrec_path)/recipes/nrutil.c

OBJS += \
./code/src/nr/gasdev.o \
./code/src/nr/gaussj.o \
./code/src/nr/spline.o \
./code/src/nr/splint.o \
./code/src/nr/nrutil.o

C_DEPS += \
./code/src/nr/gasdev.d \
./code/src/nr/gaussj.d \
./code/src/nr/spline.d \
./code/src/nr/splint.d \
./code/src/nr/nrutil.d


# Each subdirectory must supply rules for building sources it contributes
code/src/nr/%.o: $(numrec_path)/recipes/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I$(numrec_path)/recipes -I$(numrec_path)/other -O3 -c -fmessage-length=0 -std=c99 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

endif
