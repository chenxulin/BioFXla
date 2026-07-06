#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
mkdir -p results/reinvent4_fxia_pilot/raw_scale

is_finished() {
  local log="$1"
  [[ -s "${log}" ]] && grep -q "Finished REINVENT" "${log}"
}

is_running() {
  local pidfile="$1"
  [[ -s "${pidfile}" ]] && ps -p "$(cat "${pidfile}")" >/dev/null 2>&1
}

start_job() {
  local name="$1"
  local device="$2"
  local config="$3"
  local log="results/reinvent4_fxia_pilot/raw_scale/${name}.log"
  local stdout="results/reinvent4_fxia_pilot/raw_scale/${name}.stdout"
  local stderr="results/reinvent4_fxia_pilot/raw_scale/${name}.stderr"
  local pidfile="results/reinvent4_fxia_pilot/raw_scale/${name}.pid"

  if is_finished "${log}"; then
    echo "${name} finished"
    return
  fi

  if is_running "${pidfile}"; then
    echo "${name} running $(cat "${pidfile}")"
    return
  fi

  : > "${stdout}"
  : > "${stderr}"
  local reinvent_bin="${REINVENT_BIN:-reinvent}"
  local conda_lib="${CONDA_PREFIX:-/data/cxl/miniforge3/envs/aidd}/lib"
  setsid nohup env LD_LIBRARY_PATH="${conda_lib}:${LD_LIBRARY_PATH:-}" \
    "${reinvent_bin}" -d "${device}" -l "${log}" "${config}" \
    > "${stdout}" 2> "${stderr}" < /dev/null &
  echo "$!" > "${pidfile}"
  echo "${name} started $(cat "${pidfile}") ${device}"
}

start_job mol2mol_scale "${MOL2MOL_DEVICE:-cuda:0}" results/reinvent4_fxia_pilot/configs/mol2mol_scale.toml
start_job libinvent_scale "${LIBINVENT_DEVICE:-cuda:1}" results/reinvent4_fxia_pilot/configs/libinvent_scale.toml
start_job linkinvent_scale "${LINKINVENT_DEVICE:-cuda:2}" results/reinvent4_fxia_pilot/configs/linkinvent_scale.toml
start_job reinvent_denovo_scale "${REINVENT_DENOVO_DEVICE:-cuda:3}" results/reinvent4_fxia_pilot/configs/reinvent_denovo_scale.toml
