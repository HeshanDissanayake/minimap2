#!/bin/bash

INJECTOR_PATH="/home/heshds/working_dir/regws-instruction-injection"

basename=$(echo "$1" | cut -d '.' -f 1 )

dirs=("assembly" "injected" "llvm" "exec" "objdump")

for dir in "${dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        # If it doesn't exist, create the directory
        mkdir "$dir"
    fi    
done


clang  --sysroot=/opt/dev/riscv_linux_rv64g/sysroot/usr/ --gcc-toolchain=/opt/dev/riscv_linux_rv64g  -I ../zlib/zlib-1.2.11/ -DHAVE_KALLOC -O3 -S  -emit-llvm --target=riscv64 -march=rv64g  $1 -o llvm/${basename}.ll

if [ "$3" = "extra" ]; then
    /home/heshds/llvm-project/output_modified_llvm/bin/llc -O3 --march=riscv64 -mcpu=generic-rv64 -mattr=+d llvm/${basename}.ll -o assembly/${basename}_exReg.S
    python3 "$INJECTOR_PATH/inject.py" assembly/${basename}_exReg.S "$INJECTOR_PATH/encoding.json" injected/${basename}_exReg_Inj.S
    # /opt/dev/riscv_linux_rv64g/bin/riscv64-unknown-linux-gnu-gcc -static assembly/${basename}_exReg_Inj.S -o exec/${basename}_extra.o
    # /opt/dev/riscv_linux_rv64g/bin/riscv64-unknown-linux-gnu-objdump -d -S -M no-aliases exec/${basename}_extra.o > objdump/${basename}objdump_extra
else
    llc -O3 --march=riscv64 -mcpu=generic-rv64 -mattr=+d  llvm/${basename}.ll -o assembly/${basename}_normal.S
    # /opt/dev/riscv_linux_rv64g/bin/riscv64-unknown-linux-gnu-gcc -static assembly/${basename}_normal.S -o exec/${basename}_norm.o
    # /opt/dev/riscv_linux_rv64g/bin/riscv64-unknown-linux-gnu-objdump -d -S -M no-aliases exec/${basename}_norm.o > objdump/${basename}_objdump_norm
fi

