#!/bin/bash

Help() {
    echo "Setup script for genome annotation tool."
    echo
    echo "Syntax: $0 [-h|t|b|d|a]"
    echo "options:"
    echo "-h   Print this Help."
    echo "-t   Install required tools (EMBOSS, JARVIS3, GTO) and python requirements in venv."
    echo "-b   Compile C++ binaries into ./bin directory."
    echo "-a   Do both (tools + binaries)."
    echo
}

PythonVenv() { 
    python3 -m venv ${PLANT_DIR}/.venv
    source ${PLANT_DIR}/.venv/bin/activate
    pip install -r ${PLANT_DIR}/requirements.txt
}

check_and_install() {
    local pkg=$1
    local bin_check=$2
    local custom_path=$3

    if command -v "$bin_check" &> /dev/null; then
        echo "[INFO] $pkg already installed, skipping."
    else
        echo "[INFO] Installing $pkg..."
        apt-get download "$pkg" && dpkg -x ${pkg}_*.deb "$custom_path"
        echo "export PATH=$custom_path/usr/bin:\$PATH" >> ${PLANT_DIR}/env_geanno.sh
        rm -f ${pkg}_*.deb
    fi
}

Tools() {
    mkdir -p "${PLANT_DIR}/bin"
    cd "${PLANT_DIR}/bin"

    check_and_install "wget" "wget" "${PLANT_DIR}/bin/wget"
    check_and_install "make" "make" "${PLANT_DIR}/bin/make"
    check_and_install "gzip" "gunzip" "${PLANT_DIR}/bin/gzip"
    check_and_install "emboss" "getorf" "${PLANT_DIR}/bin/emboss"

    source ${PLANT_DIR}/env_geanno.sh

    # gto
    if [ ! -d "${PLANT_DIR}/bin/gto" ]; then
        git clone https://github.com/cobilab/gto.git
        (cd gto/src && make)
        echo "export PATH=${PLANT_DIR}/bin/gto/bin:\$PATH" >> ${PLANT_DIR}/env_geanno.sh
    else
        echo "[INFO] gto already installed, skipping."
    fi

    # jarvis3
    if [ ! -d "${PLANT_DIR}/bin/jarvis3" ]; then
        git clone https://github.com/cobilab/jarvis3.git
        (cd jarvis3/src && make)
        echo "export PATH=${PLANT_DIR}/bin/jarvis3/src:\$PATH" >> ${PLANT_DIR}/env_geanno.sh
    else
        echo "[INFO] jarvis3 already installed, skipping."
    fi

    cd "${PLANT_DIR}"

    echo "export PATH=${PLANT_DIR}/bin/:\$PATH" >> ${PLANT_DIR}/env_geanno.sh
    
    source ${PLANT_DIR}/env_geanno.sh

    PythonVenv;
    echo "[DONE] Tools installation complete."
}

CPP_Bin() {
    mkdir -p "${PLANT_DIR}/bin"
    g++ "${PLANT_DIR}/src/cpp/extract_characteristics_batch.cpp" -o "${PLANT_DIR}/bin/extract_characteristics_batch"
    g++ "${PLANT_DIR}/src/cpp/converter.cpp" -o "${PLANT_DIR}/bin/converter"

    echo "export PATH=${PLANT_DIR}/bin/:\$PATH" >> ${PLANT_DIR}/env_geanno.sh
    echo "[DONE] C++ binaries compiled."
}

while getopts ":htba" option; do
   case $option in
        h) Help; exit;;
        t) Tools;;
        b) CPP_Bin;;
        a) Tools; CPP_Bin;;
        \?) echo "Error: Invalid option. Use -h for help."; exit;;
   esac
done
