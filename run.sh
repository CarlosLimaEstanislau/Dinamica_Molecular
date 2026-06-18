#!/bin/bash

# ===============================
# CONFIGURAÇÃO
# ===============================
FC=gfortran
FLAGS="-O2"

SRC_DIR="src"
BUILD_DIR="build"
CONFIG_DIR="config"
PDB_DIR="pdb"
LOG_DIR="log"

MODULOS=(
  "rng_mod.f90"
  "var_mod.f90"
  "io_mod.f90"
  "utils_mod.f90"
  "forces_mod.f90"
  "ui_mod.f90"
)

PROGRAMAS=("init_conf.f90" "main.f90")
NML="$SRC_DIR/config.nml"

MASTER_DIR="$(pwd)"

# Criar diretório build se não existir
mkdir -p "$BUILD_DIR"
mkdir -p "$PDB_DIR"
mkdir -p "$CONFIG_DIR"
mkdir -p "$LOG_DIR"

altera_nml() {
    param="$1"
    valor="$2"

    if [ ! -f "$NML" ]; then
        echo "Arquivo $NML não encontrado"
        exit 1
    fi

    # 1. Normalizar whitespace (TAB -> espaço)
    sed -i.bak $'s/\t/    /g' "$NML"

    # 2. Verificar se o parâmetro existe (case-insensitive)
    if ! grep -qiE "^[[:space:]]*$param[[:space:]]*=" "$NML"; then
        echo "Parâmetro '$param' não encontrado em $NML"
        exit 1
    fi

    # 3. DETECTAR SE É STRING (tem aspas ou é texto)
    # Remove aspas simples e duplas do valor
    valor_limpo=$(echo "$valor" | sed -E "s/^['\"]//; s/['\"]$//")
    
    # Verifica se é um número (inteiro ou real)
    if [[ "$valor_limpo" =~ ^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$ ]]; then
        # É número - NÃO coloca aspas
        valor_final="$valor_limpo"
    else
        # É string - coloca aspas simples
        valor_final="'$valor_limpo'"
    fi

    # 4. Substituir valor (case-insensitive, preserva comentários)
    # Escapa caracteres especiais para o sed
    param_escaped=$(echo "$param" | sed -e 's/[\/&]/\\&/g')
    valor_escaped=$(echo "$valor_final" | sed -e 's/[\/&]/\\&/g')
    
    if sed --version >/dev/null 2>&1; then
        # GNU sed
        sed -i -E "s|(^[[:space:]]*$param_escaped[[:space:]]*=[[:space:]]*)([^!/]*)|\1$valor_escaped|I" "$NML"
    else
        # BSD sed (macOS)
        sed -i '' -E "s|(^[[:space:]]*$param_escaped[[:space:]]*=[[:space:]]*)([^!/]*)|\1$valor_escaped|I" "$NML"
    fi

    # 5. Verificar se a atualização funcionou
    if grep -qiE "^[[:space:]]*$param[[:space:]]*=[[:space:]]*$valor_final" "$NML"; then
        echo "$param atualizado para $valor_final"
    else
        echo "Falha ao atualizar $param"
        echo "Tentou colocar: $param = $valor_final"
        exit 1
    fi
}

# ===============================
# Atualizar parâmetros, se solicitado
# ===============================
if [ $# -gt 0 ]; then
    for arg in "$@"; do
        if [[ "$arg" == *=* ]]; then
            nome=${arg%%=*}
            valor=${arg#*=}
            altera_nml "$nome" "$valor"
        fi
    done
fi

# ===============================
#  Compilar módulos
# ===============================
echo -e "\e[1;33mCompilando módulos...\e[0m"

for mod in "${MODULOS[@]}"; do
    src_path="$SRC_DIR/$mod"
    
    $FC $FLAGS -c "$src_path" -J"$BUILD_DIR" -o "$BUILD_DIR/${mod%.f90}.o"
    
    if [ $? -ne 0 ]; then
        echo "Erro ao compilar módulo $mod."
        exit 1
    fi
done

# ===============================
#  Compilar e rodar programas
# ===============================
cp "$NML" "$BUILD_DIR/"

clear
echo -e "\e[1;33mSimulation starting...\e[0m"

for prog in "${PROGRAMAS[@]}"; do
    nome_exec="${prog%.f90}"
    src_path="$SRC_DIR/$prog"
    
    $FC $FLAGS "$src_path" "$BUILD_DIR"/*.o -J"$BUILD_DIR" -I"$BUILD_DIR" -o "$BUILD_DIR/$nome_exec"

    if [ $? -ne 0 ]; then
        echo "Erro ao compilar $prog."
        exit 1
    fi
    
    (cd "$BUILD_DIR" && ./"$nome_exec" "$MASTER_DIR")
done

echo -e "\e[1;33mEnd of simulation!\e[0m"