#!/bin/bash
set -e

DATA_DIR="/data"
GM_KEY="/home/hacga_user/.gm_key"

if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: $DATA_DIR is not mounted or does not exist."
    exit 1
fi

if [ ! -f "$GM_KEY" ]; then
    echo "ERROR: $GM_KEY is not mounted!"
    echo "       Please mount your .gm_key file to $GM_KEY in the container."
    exit 1
fi

HOST_UID=$(stat -c "%u" "$DATA_DIR")
HOST_GID=$(stat -c "%g" "$DATA_DIR")

if [ "$HOST_UID" -eq 0 ]; then
    echo "WARNING: $DATA_DIR is owned by root. Using default UID 1000."
    HOST_UID=1000
fi
if [ "$HOST_GID" -eq 0 ]; then
    HOST_GID=1000
fi

echo "Using UID=$HOST_UID, GID=$HOST_GID"

if ! getent group "$HOST_GID" >/dev/null 2>&1; then
    groupadd -g "$HOST_GID" hacga_group
fi

if ! id -u hacga_user >/dev/null 2>&1; then
    useradd -o -m -u "$HOST_UID" -g "$HOST_GID" -s /bin/bash hacga_user
else
    usermod -o -u "$HOST_UID" -g "$HOST_GID" hacga_user 2>/dev/null || true
fi

if [ ! -r "$GM_KEY" ]; then
    echo "ERROR: $GM_KEY is not readable!"
    exit 1
fi

chown hacga_user /home/hacga_user
chown hacga_user "$GM_KEY"
chmod 600 "$GM_KEY"

chown hacga_user "$DATA_DIR"

echo "=== HACGA1 Container ==="
echo "Running as 'hacga_user' (UID=$HOST_UID, GID=$HOST_GID)"

source /opt/conda/etc/profile.d/conda.sh

if [ $# -eq 0 ]; then
    exec gosu hacga_user bash
else
    exec gosu hacga_user conda run -n ann --no-capture-output "$@"
fi
